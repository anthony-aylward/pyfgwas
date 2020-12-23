#===============================================================================
# pyfgwas.py
#===============================================================================

"""Apply the workflow suggested by the fgwas manual"""



# Imports ======================================================================

import argparse
import copy
import gzip
import itertools
import os
import operator
import socket
import subprocess
import sys
import tempfile
import time

from multiprocessing import Pool

from pyfgwas.fgwasplot import plot_params




# Constants ====================================================================

PATH = os.environ.get('PYFGWAS_PATH', '/usr/local/bin/fgwas')




# Class definitions ============================================================

class InputFileHeader():
    """The header of the input file
    """

    non_annotation_column_names = {
            'SNPID',
            'CHR',
            'POS',
            'Z',
            'F',
            'N',
            'NCASE',
            'NCONTROL',
            'LNBF',
            'SEGNUMBER'
        }
    
    def __init__(self, args):
        self.args = args
        with gzip.open(args.input, 'rt') as f:
            self.tuple = tuple(f.readline().replace('\n', '').split(' '))
        self.annotations = {
            column_name
            for column_name in self.tuple
            if column_name not in self.non_annotation_column_names
        }
        self.validate()
    
    def validate(self):
        """Validate the input file header
        """
        
        for column_name in self.non_annotation_column_names:
            if self.tuple.count(column_name) > 1:
                raise SystemExit(
                    (
                        'ERROR: Invalid header on input file: Either {0} is '
                        'duplicated, or there is an annotation called {0}.'
                    ).format(column_name)
                )
            elif (
                (
                    column_name
                    not in {'N', 'NCASE', 'NCONTROL', 'LNBF', 'SEGNUMBER'}
                )
                and (column_name not in self.tuple)
            ):
                raise SystemExit(
                    'ERROR: Invalid header on input file: {} is missing.'
                    .format(column_name)
                )
        if {'N', 'NCASE', 'NCONTROL'} <= set(self.tuple):
            raise SystemExit(
                'ERROR: N, NCASE, and NCONTROL were all found in the input '
                'file header. Please use EITHER quantitative OR case-control '
                'format.'
            )


class IndividualAnnotationResults():
    """Individual annotation results
    """
    
    def __init__(self, args, header):
        self.args = args
        self.annotations = header.annotations
        self.header = header
        self.data = None
        self.defined_ci_annotations = set()
    
    def collect(self, base_model):
        """Collect individual annotation results
        """
        
        with Pool(processes=min(self.args.processes, len(self.annotations))) as pool:
            self.data = (
                (
                    (
                        'parameter',
                        'ln(lk)',
                        'CI_lo',
                        'estimate',
                        'CI_hi',
                        'converged'
                    ),
                )
                + tuple(
                    sorted(
                        pool.starmap(
                            collect_individual_annotation_result,
                            (
                                (self.args, annotation, self.header, base_model)
                                for annotation in (
                                    self.annotations
                                    - set(self.args.base_model.split('+'))
                                )
                            )
                        ),
                        key=operator.itemgetter(1),
                        reverse=True
                    )
                )
            )
    
    def export(self, output_path):
        """Export individual annotation results"""
        
        if self.data:
             with open(output_path, 'w') as f:
                f.write(
                    '\n'.join(
                        '\t'.join((str(entry) for entry in line))
                        for line in self.data
                    )
                    + '\n'
                )
        else:
            raise RuntimeError(
                'Cannot export individual annotation results with no data '
                'present.'
            )
    
    def identify_defined_ci_annotations(self):
        if self.data:
            for parameter, llk, ci_lo, estimate, ci_hi, converged in (
                self.data[1:]
            ):
                try:
                    float(ci_lo)
                    float(ci_hi)
                    if (
                        ((float(ci_lo) > 0) or (float(ci_hi) < 0))
                        if self.args.strict_joint_model
                        else True
                    ):
                        self.defined_ci_annotations.add(parameter)
                except ValueError:
                    pass
        if not self.defined_ci_annotations:
            raise SystemExit(
                'None of the provided annotations were viable (all '
                'individual results had undefined confidence intervals). '
                'Terminating analysis.'
            )
    
    def identify_best_starting_state(self):
        if self.data:
            best_individual_llk = max(
                float(llk)
                for parameter, llk, ci_lo, estimate, ci_hi, converged
                in self.data[1:]
            )
            best_individual_annotations = []
            for parameter, llk, ci_lo, estimate, ci_hi, converged in (
                self.data[1:]
            ):
                if (
                    float(llk) == best_individual_llk and (
                        (
                            parameter.replace('_ln', '')
                            in self.defined_ci_annotations
                        ) if self.defined_ci_annotations else True
                    )
                ):
                    best_individual_annotations.append(
                        (parameter, llk, ci_lo, estimate, ci_hi)
                    )
            if len(best_individual_annotations) == 1:
                best_individual_annotation = best_individual_annotations[0]
                estimates = collect_estimates(
                    '{}-indiv.tsv'.format(self.args.output),
                    self.args.base_model.split('+')
                    + [best_individual_annotation[0]]
                )
                self.best_starting_state = {
                    'annotations': (
                        bool(self.args.base_model)
                        * self.args.base_model.split('+')
                        + [best_individual_annotation[0]]
                    ),
                    'llk': best_individual_annotation[1],
                    'estimates': estimates,
                    'xvl': float('-inf'),
                    'xv_penalty': 0,
                    'output_files': {}
                }
            elif len(best_individual_annotations) > 1:
                annotation_combinations = get_annotation_combinations(
                    (
                        parameter
                        for parameter, llk, ci_lo, estimate, ci_hi
                        in best_individual_annotations
                    )
                )
                with tempfile.TemporaryDirectory() as temp_dir_name:
                    with Pool(
                        processes=min(
                            self.args.processes,
                            len(annotation_combinations)
                        )
                    ) as pool:
                        pool.map(
                            call_fgwas,
                            (
                                {
                                    'args': self.args,
                                    'annotation': '+'.join(
                                        (self.args.base_model,)
                                        + tuple(annotation_combination)
                                    ),
                                    'tag': '+'.join(annotation_combination),
                                    'header': self.header,
                                    'output_dir': temp_dir_name
                                }
                                for annotation_combination
                                in annotation_combinations
                            )
                        )
                    (
                        best_combo_llk,
                        best_combos
                    ) = best_annotation_combinations_llk(
                        annotation_combinations,
                        temp_dir_name
                    )
                    if (
                        (best_combo_llk > best_individual_llk)
                        and len(best_combos) == 1
                    ):
                        best_combo = list(best_combos[0])
                        estimates = collect_estimates(
                            '{}/{}.params'
                            .format(temp_dir_name, '+'.join(best_combo)),
                            self.annotations
                        )
                        self.best_starting_state = {
                            'annotations': (
                                bool(self.args.base_model)
                                * self.args.base_model.split('+')
                                + best_combo
                            ),
                            'llk': best_combo_llk,
                            'estimates': estimates,
                            'xvl': float('-inf'),
                            'xv_penalty': 0,
                            'output_files': {}
                        }
                    else:
                        print(
                            'Next best annotation is ambiguous, taking a null '
                            'step and proceeding to cross-validation phase'
                        )
        else:
            raise Exception(
                'Cannot identify annotations with with well-defined confidence '
                'intervals with no data present.'
            )
        

class FgwasModel():
    """An FGWAS model
    """
    
    def __init__(self, args, header):
        self.args = args
        self.header = header
        self.annotations = []
        self.llk = float('-inf')
        self.estimates = {}
        self.xvl = float('-inf')
        self.cross_validation_penalty = 0
        self.output_files = {}
        if args.base_model:
            self.annotations = args.base_model.split('+')
            with tempfile.TemporaryDirectory() as temp_dir_name:
                call_fgwas(
                    {
                        'args': args,
                        'annotation': args.base_model,
                        'header': header,
                        'output_dir': temp_dir_name,
                        'tag': args.base_model
                    }
                )
                self.update(self.args.base_model, temp_dir_name)
        if args.resume:
            self.annotations.extend(args.resume.split('+'))
            with tempfile.TemporaryDirectory() as temp_dir_name:
                call_fgwas(
                    {
                        'args': args,
                        'annotation': '+'.join(self.annotations),
                        'header': header,
                        'output_dir': temp_dir_name,
                        'tag': '+'.join(self.annotations)
                    }
                )
                self.update('+'.join(self.annotations), temp_dir_name)
        self.cache()
    
    def update(self, tag, temp_dir_name):
        with open('{}/{}.llk'.format(temp_dir_name, tag), 'r') as f:
            self.llk = float(f.readline().rstrip('\n').replace('ln(lk): ', ''))
        with open('{}/{}.params'.format(temp_dir_name, tag), 'r') as f:
            for line in f:
                parsed_line = tuple(line.rstrip('\n').split(' '))
                parameter = parsed_line[0][:-3]
                for annotation in self.annotations:
                    if annotation == parameter:
                        self.estimates[annotation] = tuple(parsed_line[1:])
                        break
        with open('{}/{}.stderr'.format(temp_dir_name, tag), 'r') as f:
            self.converged = not (
                f.readline().rstrip('\n') == 'WARNING: failed to converge'
            )
    
    def clear(self):
        """Clear the model to an empty state
        """
        
        self.cache()
        self.annotations = []
        self.llk = float('-inf')
        self.estimates = {}
        self.xvl = float('-inf')
        self.cross_validation_penalty = 0
        self.output_files = {}
    
    def cache(self):
        """Cache the model state
        """
        
        self.annotations_cache = list(self.annotations)
        self.llk_cache = float(self.llk)
        self.estimates_cache = dict(self.estimates)
        self.xvl_cache = float(self.xvl)
        self.cross_validation_penalty_cache = float(
            self.cross_validation_penalty
        )
        self.output_files_cache = dict(self.output_files)
    
    def revert(self):
        """Revert to the cached model state
        """
        
        self.annotations = self.annotations_cache
        self.llk = self.llk_cache
        self.estimates = self.estimates_cache
        self.xvl = self.xvl_cache
        self.cross_validation_penalty = self.cross_validation_penalty_cache
        self.output_files = self.output_files_cache
        print('Reverted last step')
    
    def set_state(self, state):
        """Set the model state according to a dictionary definition
        """
        
        self.cache()
        self.annotations = state['annotations']
        self.llk = state['llk']
        self.estimates = state['estimates']
        self.xvl = state['xvl']
        self.cross_validation_penalty = state['xv_penalty']
        self.output_files = state['output_files']
    
    def collect_output_files(self, annotation, temp_dir_name):
        '''
        Collect the output files tagged with a specific annotation
        '''
        for extension in ('llk', 'params', 'ridgeparams', 'stderr'):
            with open(
                '{}/{}.{}'.format(temp_dir_name, annotation, extension),
                'r'
            ) as f:
                self.output_files[extension] = f.read()
        for extension in ('bfs', 'segbfs'):
            with gzip.open(
                '{}/{}.{}.gz'.format(temp_dir_name, annotation, extension),
                'rb'
            ) as f:
                self.output_files[extension] = f.read()
    
    def export(self, tag):
        """Export the output file state
        """
        
        for extension in ('llk', 'params', 'ridgeparams', 'stderr'):
            with open(
                '{}-{}.{}'.format(self.args.output, tag, extension),
                'w'
            ) as f:
                f.write(self.output_files[extension])
        for extension in ('bfs', 'segbfs'):
            with gzip.open(
                '{}-{}.{}.gz'.format(self.args.output, tag, extension),
                'wb'
            ) as f:
                f.write(self.output_files[extension])
        if self.args.generate_plots:
            plot_params(
                '{}-{}.params'.format(self.args.output, tag),
                '{}-{}.pdf'.format(self.args.output, tag),
                title=self.args.plot_title.format(tag)
            )
    
    def report_appended_annotation(self, annotation, llk):
        print(
            'Added {} to joint model (llk: {})'
            .format(annotation, str(llk))
        )
        with open('{}.cache'.format(self.args.output), 'w') as f:
            f.write('+'.join(self.annotations))
    
    def append(self, annotation):
        """Append an annotation to the model and determine the resulting joint
        model
        """
        
        self.cache()
        self.annotations.append(annotation)
        with tempfile.TemporaryDirectory() as temp_dir_name:
            call_fgwas(
                {
                    'args': self.args,
                    'annotation': '+'.join(self.annotations),
                    'header': self.header,
                    'output_dir': temp_dir_name
                }
            )
            self.update(annotation, temp_dir_name)
        
    def append_best_annotation(self, individual_results):
        """Identify the next best annotation, append it to the model, and
        determine the resulting joint model
        """
        
        self.cache()
        remaining_annotations = {
            annotation
            for annotation in individual_results.defined_ci_annotations
            if annotation not in self.annotations
        }
        with tempfile.TemporaryDirectory() as temp_dir_name:
            with Pool(
                processes=min(
                    self.args.processes,
                    max(len(remaining_annotations), 1),
                )
            ) as pool:
                pool.map(
                    call_fgwas,
                    (
                        {
                            'args': self.args,
                            'annotation': '+'.join(
                                self.annotations + [annotation]
                            ),
                            'header': individual_results.header,
                            'output_dir': temp_dir_name
                        }
                        for annotation in remaining_annotations
                    )
                )
            llk_list = []
            individual_estimates = {
                parameter: (
                    float(
                        '-inf'
                        if ci_lo.startswith('<')
                        else ci_lo.replace('fail', '-inf')
                    ),
                    float(ci_hi.replace('fail', 'inf'))
                )
                for parameter, llk, ci_lo, estimate, ci_hi, converged
                in individual_results.data[1:]
                if parameter
                in individual_results.defined_ci_annotations
            }
            for annotation in remaining_annotations:
                consistent_estimates = True
                with open(
                    '{}/{}.params'.format(temp_dir_name, annotation)
                ) as f:
                    f.readline()
                    for line in f:
                        parameter, ci_lo, estimate, ci_hi = (
                            line.rstrip('\n').split(' ')
                        )
                        if parameter.replace('_ln', '') != 'pi_region':
                            if (
                                (
                                    float(
                                        '-inf'
                                        if ci_lo.startswith('<')
                                        else ci_lo.replace('fail', '-inf')
                                    )
                                    > individual_estimates[
                                        parameter.replace('_ln', '')
                                    ][1]
                                    or
                                    float(ci_hi.replace('fail', 'inf'))
                                    < individual_estimates[
                                        parameter.replace('_ln', '')
                                    ][0]
                                )
                                and (
                                    parameter[:-3]
                                    not in self.args.base_model.split('+')
                                )
                            ):
                                consistent_estimates = False
                if (
                    consistent_estimates
                    if individual_results.args.consistent_estimates_only
                    else True
                ):
                    with open(
                        '{}/{}.llk'.format(temp_dir_name, annotation)
                    ) as f0, open(
                        '{}/{}.stderr'.format(temp_dir_name, annotation)
                    ) as f1:
                        llk_list.append(
                            (
                                float(
                                    f0
                                    .readline()
                                    .replace('ln(lk): ', '')
                                    .rstrip('\n')
                                ),
                                annotation,
                                (
                                    f1.readline().rstrip('\n')
                                    != 'WARNING: failed to converge'
                                )
                            )
                        )
            llk_list = [
                (llk, annotation)
                for llk, annotation, converged in llk_list
#                 if converged
            ]
            if individual_results.args.consistent_estimates_only:
                print(
                    '{} annotations can be added with consistent parameter '
                    'estimates'
                    .format(len(llk_list))
                )
            if len(llk_list) == 0:
                print('No more annotations can be added')
            else:
                best_llk = max(llk for llk, annotation in llk_list)
                best_annotations = tuple(
                    annotation
                    for llk, annotation in llk_list
                    if llk == best_llk
                )
                if len(best_annotations) == 1:
                    best_annotation = best_annotations[0]
                    self.llk = best_llk
                    self.annotations.append(best_annotation)
                    with open(
                        '{}/{}.params'.format(temp_dir_name, best_annotation),
                        'r'
                    ) as f:
                        for line in f:
                            parsed_line = tuple(
                                line.replace('\n', '').split(' ')
                            )
                            parameter = parsed_line[0][:-3]
                            for annotation in self.annotations:
                                if annotation == parameter:
                                    self.estimates[annotation] = tuple(
                                        parsed_line[1:]
                                    )
                                    break
                    self.collect_output_files(best_annotation, temp_dir_name)
                    self.report_appended_annotation(best_annotation, best_llk)
                elif len(best_annotations) > 1:
                    annotation_combinations = get_annotation_combinations(
                        best_annotations
                    )
                    with Pool(
                        processes=min(
                            self.args.processes,
                            len(annotation_combinations)
                        )
                    ) as pool:
                        pool.map(
                            call_fgwas,
                            (
                                {
                                    'args': self.args,
                                    'annotation':
                                        '+'.join(
                                            self.annotations
                                            + list(annotation_combination)
                                        ),
                                    'tag': '+'.join(annotation_combination),
                                    'header': individual_results.header,
                                    'output_dir': temp_dir_name
                                }
                                for annotation_combination
                                in annotation_combinations
                            )
                        )
                    best_combo_llk, best_combos = best_annotation_combinations_llk(
                        annotation_combinations,
                        temp_dir_name
                    )
                    if (best_combo_llk > best_llk) and len(best_combos) == 1:
                        best_combo = list(best_combos[0])
                        self.llk = best_combo_llk
                        self.annotations.extend(best_combo)
                        with open(
                            '{}/{}.params'
                            .format(temp_dir_name, '+'.join(best_combo))
                        ) as f:
                            for line in f:
                                parsed_line = tuple(
                                    line.replace('\n', '').split(' ')
                                )
                                parameter = parsed_line[0][:-3]
                                for annotation in self.annotations:
                                    if annotation == parameter:
                                        self.estimates[annotation] = tuple(
                                            parsed_line[1:]
                                        )
                                        break
                        self.collect_output_files(
                            '+'.join(best_combo),
                            temp_dir_name
                        )
                        self.report_appended_annotation(
                            '+'.join(best_combo),
                            best_llk
                        )
                    else:
                        print(
                            'Next best annotation is ambiguous, taking a null '
                            'step and proceeding to cross-validation phase'
                        )

    def calibrate_cross_validation_penalty(self, header, processes=1):
        """Calibrate the cross validation penalty
        """
        
        xv_penalties = tuple(
            0.5 * x / max(processes, 10)
            for x in range(1, max(processes, 10) + 1)
        )
        with tempfile.TemporaryDirectory() as temp_dir_name:
            with Pool(processes=processes) as pool:
                pool.map(
                    call_fgwas,
                    tuple(
                        {
                            'args': self.args,
                            'annotation': '+'.join(self.annotations),
                            'header': header,
                            'output_dir': temp_dir_name,
                            'xv_penalty': xv_penalty
                        }
                        for xv_penalty in xv_penalties
                    )
                )
            xv_likelihoods = []
            for xv_penalty in xv_penalties:
                with open(
                    '{}/p{}.ridgeparams'.format(temp_dir_name, str(xv_penalty)),
                    'r'
                ) as f:
                    xv_likelihoods.append(
                        (
                            float(
                                f
                                .read()
                                .split('\n')[-2]
                                .replace('X-validation penalize ln(lk): ', '')
                                .replace('\n', '')
                            ),
                            xv_penalty
                        )
                    )
        best_xvl_penalty = (
            sorted(xv_likelihoods, key=operator.itemgetter(0), reverse=True)
            [0]
        )
        self.xvl = best_xvl_penalty[0]
        self.cross_validation_penalty = best_xvl_penalty[1]
        
    def remove_worst_annotation(self, header):
        """Identify the annotation contributing least to the model and remove it
        """
        
        if len(self.annotations) < 2:
            raise Exception(
                'Can\'t remove annotations from a model with only one '
                'annotation.'
            )
        self.cache()
        with tempfile.TemporaryDirectory() as temp_dir_name:
            with Pool(
                processes=min(self.args.processes,len(self.annotations))
            ) as pool:
                pool.map(
                    call_fgwas,
                    (
                        {
                            'args': self.args,
                            'annotation': annotation_0,
                            'header': header,
                            'output_dir': temp_dir_name,
                            'xv_penalty': self.cross_validation_penalty,
                            'drop': '+'.join(
                                annotation_1
                                for annotation_1 in self.annotations
                                if (annotation_0 != annotation_1)
                            )
                        }
                        for annotation_0 in self.annotations
                    )
                )
            xvl_list = []
            for annotation in self.annotations:
                with open(
                    '{}/{}.ridgeparams'.format(temp_dir_name, annotation),
                    'r'
                ) as f:
                    xvl_list.append(
                        (
                            float(
                                f
                                .read()
                                .split('\n')[-2]
                                .replace('X-validation penalize ln(lk): ', '')
                                .replace('\n', '')
                            ),
                            annotation
                        )
                    )
            best_xvl = max(xvl for xvl, annotation in xvl_list)
            worst_annotations = tuple(
                annotation
                for xvl, annotation in xvl_list
                if xvl == best_xvl
            )
            if len(worst_annotations) == 1:
                worst_annotation = worst_annotations[0]
                self.xvl = best_xvl
                self.annotations.remove(worst_annotation)
                with open(
                    '{}/{}.params'.format(temp_dir_name, worst_annotation)
                ) as f:
                    for line in f:
                        parsed_line = tuple(line.replace('\n', '').split(' '))
                        parameter = parsed_line[0][:-3]
                        for annotation in self.annotations:
                            if annotation == parameter:
                                self.estimates[annotation] = tuple(
                                    parsed_line[1:]
                                )
                                break
                self.collect_output_files(worst_annotation, temp_dir_name)
                print(
                    'Dropped {} from joint model (xvl: {})'
                    .format(worst_annotation, str(best_xvl))
                )
            elif len(worst_annotations) > 1:
                annotation_combinations = get_annotation_combinations(
                    worst_annotations
                )
                with Pool(
                    processes=min(
                        self.args.processes,
                        len(annotation_combinations)
                    )
                ) as pool:
                    pool.map(
                        call_fgwas,
                        (
                            {
                                'args': self.args,
                                'annotation': '+'.join(
                                    self.annotations
                                    + list(annotation_combination)
                                ),
                                'tag': '+'.join(sorted(annotation_combination)),
                                'header': header,
                                'output_dir': temp_dir_name,
                                'xv_penalty': self.cross_validation_penalty,
                                'drop': '+'.join(
                                    annotation
                                    for annotation in self.annotations
                                    if
                                    (annotation not in annotation_combination)
                                )
                            }
                            for annotation_combination
                            in annotation_combinations
                        )
                    )
                best_combo_xvl, best_combos = best_annotation_combinations_xvl(
                    annotation_combinations,
                    temp_dir_name
                )
                if best_combo_xvl >= best_xvl:
                    best_combo = sorted(
                        list(
                            {
                                annotation
                                for combo in best_combos
                                for annotation in combo
                            }
                        )
                    )
                    self.xvl = best_combo_xvl
                    for annotation in best_combo:
                        self.annotations.remove(annotation)
                    with open(
                        '{}/{}.params'
                        .format(temp_dir_name, '+'.join(best_combo))
                    ) as f:
                        for line in f:
                            parsed_line = tuple(
                                line.replace('\n', '').split(' ')
                            )
                            parameter = parsed_line[0][:-3]
                            for annotation in self.annotations:
                                if annotation == parameter:
                                    self.estimates[annotation] = tuple(
                                        parsed_line[1:]
                                    )
                                    break
                    self.collect_output_files(
                        '+'.join(best_combo),
                        temp_dir_name
                    )
                    print(
                        'Dropped {} from joint model (xvl: {})'
                        .format('+'.join(best_combo), str(best_xvl))
                    )
                else:
                    print(
                        'Next annotation to drop is ambiguous, taking a null '
                        'step and terminating analysis'
                    )
        
        
        

# Functions ====================================================================

def call_fgwas(args_dict):
    """A function for calling FGWAS.
    """
    args = args_dict['args']
    annotation = args_dict['annotation']
    header = args_dict['header']
    output_dir = args_dict['output_dir']
    tag = (
        args_dict['tag']
        if 'tag' in args_dict.keys()
        else annotation.split('+')[-1]
    )
    xv_penalty = (
        args_dict['xv_penalty']
        if 'xv_penalty' in args_dict.keys()
        else None
    )
    drop = args_dict['drop'] if 'drop' in args_dict.keys() else ''
    fgwas_command_line = (
        (
            'fgwas', '-print',
            '-i', args.input,
            '-w', drop if bool(drop) else annotation,
            '-o', (
                '{}/{}'.format(output_dir, tag)
                if (not bool(xv_penalty)) or bool(drop)
                else '{}/p{}'.format(output_dir, str(xv_penalty))
            )
        )
        + bool(xv_penalty) * ('-xv', '-p', str(xv_penalty))
        + bool(args.bed) * ('-bed', args.bed)
        + ({'NCASE', 'NCONTROL'} <= set(header.tuple)) * ('-cc',)
        + ('SEGNUMBER' in header.tuple) * ('-fine',)
        + bool(args.window) * ('-k', args.window)
    )
    with open('{}/{}.stderr'.format(output_dir, tag), 'w') as f:
        subprocess.call(
            fgwas_command_line,
            stdout=subprocess.DEVNULL,
            stderr=f
        )


def collect_individual_annotation_result(args, annotation, header, base_model):
    """Collect the result for an individual annotation
    """
    
    model = copy.deepcopy(base_model)
    if args.resume and not args.base_model:
        model.clear()
    elif args.resume and args.base_model:
        temp_args = args
        temp_args.resume = None
        model = FgwasModel(temp_args, header)
    model.append(annotation)
    return (
        (annotation, model.llk)
        + model.estimates[annotation]
        + (str(model.converged),)
    )


def collect_estimates(filename, annotations):
    estimates = {}
    with open(filename) as f:
        for line in f:
            parsed_line = tuple(
                line.rstrip('\n').split(' ')
            )
            parameter = parsed_line[0][:-3]
            for annotation in annotations:
                if annotation == parameter:
                    estimates[annotation] = tuple(
                        parsed_line[1:]
                    )
                    break
    return estimates


def get_annotation_combinations(annotation_sequence):
    return tuple(
        itertools.chain.from_iterable(
            itertools.combinations(annotation_sequence, length)
            for length in range(2, len(annotation_sequence) + 1)
        )
    )


def best_annotation_combinations_llk(annotation_combinations, dir_name):
    llk_list = []
    for annotation_combination in annotation_combinations:
        with open(
            '{}/{}.llk'.format(dir_name, '+'.join(annotation_combination))
        ) as f:
            llk_list.append(
                (
                    float(
                        f
                        .readline()
                        .replace('ln(lk): ', '')
                        .replace('\n', '')
                    ),
                    annotation_combination
                )
            )
    best_combo_llk = max(llk for llk, annotation_combination in llk_list)
    best_combos = tuple(
        annotation_combination
        for llk, annotation_combination in llk_list
        if llk == best_combo_llk
    )
    return best_combo_llk, best_combos


def best_annotation_combinations_xvl(annotation_combinations, dir_name):
    xvl_list = []
    for annotation_combination in annotation_combinations:
        with open(
            '{}/{}.ridgeparams'
            .format(dir_name, '+'.join(sorted(annotation_combination)))
        ) as f:
            xvl_list.append(
                (
                    float(
                        f
                        .read()
                        .split('\n')[-2]
                        .replace('X-validation penalize ln(lk): ', '')
                        .replace('\n', '')
                    ),
                    annotation_combination
                )
            )
    best_combo_xvl = max(
        xvl for xvl, annotation_combination in xvl_list
    )
    best_combos = [
        annotation_combination
        for xvl, annotation_combination in xvl_list
        if xvl == best_combo_xvl
    ]
    return best_combo_xvl, best_combos


def main():
    args = parse_arguments()
    assert os.path.isfile(PATH), (
        'fgwas not found, please set environment variable PYFGWAS_PATH'
    )
    print('Reading input file header')
    header = InputFileHeader(args)
    print('Initializing model')
    model = FgwasModel(args, header)
    individual_results = IndividualAnnotationResults(args, header)
    print('Collecting individual annotation results')
    individual_results.collect(model)
    print('Exporting individual annotation results')
    individual_results.export(f'{args.output}-indiv.tsv')
    print('Identifying annotations with well-defined confidence intervals')
    individual_results.identify_defined_ci_annotations()
    print(
        '{} annotations with well-defined confidence intervals.'
        .format(len(individual_results.defined_ci_annotations))
    )
    if not args.resume:
        individual_results.identify_best_starting_state()
        model.set_state(individual_results.best_starting_state)
    print(
        'Constructing joint model, beginning with: {} (llk: {})'
        .format(', '.join(model.annotations), model.llk)
    )
    while (
        (
            individual_results
            .defined_ci_annotations
            .union(set(args.base_model.split('+')))
        )
        > set(model.annotations)
    ):
        model.append_best_annotation(individual_results)
        if model.llk <= model.llk_cache + args.threshold:
            model.revert()
            break
    print('Exporting pre-cross-validation results')
    model.export('pre-xv')
    model.calibrate_cross_validation_penalty(header, processes=args.processes)
    print('Beginning cross-validation phase')
    number_of_annotations = len(model.annotations)
    for iteration in range(number_of_annotations - 1):
        model.remove_worst_annotation(header)
        if model.xvl <= model.xvl_cache:
            model.revert()
            break
    print('Exporting post-cross-validation results')
    model.export('post-xv')
    print('Workflow complete.')


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            'Apply the workflow suggested by the fgwas manual'
        )
    )
    parser.add_argument(
        'input',
        metavar='<path/to/fgwas/input_file.txt.gz>',
        help='Path to fgwas input file'
    )
    parser.add_argument(
        'output',
        metavar='<prefix/for/output_files>',
        help='Prefix for output files'
    )
    resource_group = parser.add_argument_group('resource arguments')
    resource_group.add_argument(
        '-p',
        '--processes',
        metavar='<int>',
        default=1,
        type=int,
        help='Maximum number of fgwas processes to launch'
    )
    parameter_group = parser.add_argument_group('parameter arguments')
    parameter_group.add_argument(
        '-b',
        '--bed',
        metavar='<path/to/bed_file.bed>',
        help='Path to .bed file containing region definitions'
    )
    parameter_group.add_argument(
        '-k',
        '--window',
        metavar='<int>',
        help='Window size'
    )
    workflow_group = parser.add_argument_group('workflow arguments')
    workflow_group.add_argument(
        '-t',
        '--threshold',
        metavar='<float>',
        type=float,
        default=0,
        help='Likelihood threshold for model improvement'
    )
    workflow_group.add_argument(
        '--base-model',
        default='',
        metavar='<list+of+annotations>',
        help='"+"-separated list of annotations to use as a base model.'
    )
    workflow_group.add_argument(
        '--resume',
        default='',
        metavar='<list+of+annotations>',
        help=(
            '"+"-separated list of annotations to include, for resuming '
            'generation of a joint model.'
        )
    )
    alternate_workflow_group = parser.add_argument_group(
        'alternate workflow arguments'
    )
    alternate_workflow_group.add_argument(
        '--consistent-estimates-only',
        action='store_true',
        help=(
            'Add an annotation to the joint model only if all parameter '
            'estimates are consistent with their individual estimates.'
        )
    )
    alternate_workflow_group.add_argument(
        '--strict-joint-model',
        action='store_true',
        help=(
            'When building a joint model, use only annotations whose '
            'individual confidence intervals do not overlap zero.'
        )
    )
    alternate_workflow_group.add_argument(
        '--generate-plots',
        action='store_true',
        help='generate plots for joint models'
    )
    alternate_workflow_group.add_argument(
        '--plot-title',
        metavar='<"plot title">',
        default='',
        help=(
            'Provide a title for the plots generated when --generate-plots is '
            'on'
        )
    )
    args = parser.parse_args()
    assert args.threshold >= 0, 'The likelihood threshold must be nonnegative.'
    return args
