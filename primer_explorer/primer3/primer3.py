import subprocess
from primer_explorer import config
from itertools import zip_longest


def reverse_complement(seq):
    return ''.join([config.COMPLEMENTARY_BASES.get(nucl, 'N')
                    for nucl in reversed(seq.decode())]).encode()


def _get_parameters_for_pair():
    parameters = "SEQUENCE_PRIMER_REVCOMP={8}\n"
    parameters += "PRIMER_PAIR_MAX_COMPL_ANY_TH={5}\n"
    parameters += "PRIMER_PAIR_MAX_COMPL_END_TH={5}\n"
    parameters += "PRIMER_PAIR_MAX_COMPL_END={6}\n"
    parameters += "PRIMER_PAIR_MAX_COMPL_ANY={7}\n"
    # 40 by default, you should use your required temp
    parameters += "PRIMER_PAIR_MAX_DIFF_TM={9}\n"
    return parameters


def _get_parameters_for_kmer():
    parameters = "SEQUENCE_ID=Check primers\n"
    parameters += "PRIMER_TASK=check_primers\n"
    parameters += "PRIMER_THERMODYNAMIC_PARAMETERS_PATH={0}\n"
    parameters += "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT={1}\n"
    parameters += "PRIMER_MAX_NS_ACCEPTED=3\n"
    parameters += "PRIMER_NUM_RETURN=1\n"
    parameters += "PRIMER_EXPLAIN_FLAG=1\n"
    parameters += "PRIMER_MIN_TM={2}\n"
    parameters += "PRIMER_MIN_SIZE={3}\n"
    parameters += "PRIMER_PICK_RIGHT_PRIMER=0\n"
    parameters += "PRIMER_PICK_INTERNAL_OLIGO=0\n"
    parameters += "PRIMER_PICK_LEFT_PRIMER=0\n"
    parameters += "PRIMER_NUM_RETURN=5\n"
    parameters += "SEQUENCE_PRIMER={4}\n"
    parameters += "PRIMER_PICK_ANYWAY=0\n"
    parameters += "PRIMER_MAX_SELF_ANY_TH={5}\n"
    parameters += "PRIMER_MAX_SELF_END_TH={5}\n"
    parameters += "PRIMER_MAX_SELF_END={6}\n"
    parameters += "PRIMER_MAX_SELF_ANY={7}\n"
    parameters += "PRIMER_MAX_HAIRPIN_TH={6}\n"
    parameters += "PRIMER_INTERNAL_MAX_HAIRPIN_TH={6}\n"
    return parameters


def _run_primer3(primer3_binary, arguments):
    process = subprocess.run(primer3_binary, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             input=arguments.encode('utf-8'))

    stdout = process.stdout
    stderr = process.stderr

    if process.returncode == 0:
        results = stdout.decode("utf-8").split("=\n")[:-1]
        return results

    else:
        msg = "Program error: {}".format(stderr)
        raise RuntimeError("{}".format(msg))


def _run_primer3_with_kmer(kmer, filtering_criteria):
    parameters = _get_parameters_for_kmer()
    parameters += "=\n"
    arguments = ""
    max_comp_th = filtering_criteria["tm_threshold"] - 10

    if max_comp_th < 0:
        max_comp_th = 0

    for thermodynamic_mode in (0, 1):
        min_len = len(kmer)
        arguments += parameters.format(filtering_criteria["primer3_config_fpath"],
                                       thermodynamic_mode,
                                       filtering_criteria["tm_threshold"],
                                       min_len, kmer, max_comp_th,
                                       filtering_criteria["max_complementary_end"],
                                       filtering_criteria["max_complementary_any"])
    return _run_primer3(filtering_criteria["primer3_core_fpath"], arguments)


def _run_primer3_with_pair(pair, filtering_criteria):
    parameters = _get_parameters_for_kmer() + _get_parameters_for_pair()
    parameters += "=\n"

    arguments = ""
    max_comp_th = filtering_criteria["tm_threshold"] - 10

    if max_comp_th < 0:
        max_comp_th = 0
    fwd_primer = pair[0].decode()
    rev_primer = reverse_complement(pair[1].encode()).decode()
    min_len = len(fwd_primer)
    for thermodynamic_mode in (0, 1):
        arguments += parameters.format(filtering_criteria["primer3_config_fpath"],
                                       thermodynamic_mode,
                                       filtering_criteria["tm_threshold"],
                                       min_len, fwd_primer, max_comp_th,
                                       filtering_criteria["max_complementary_end"],
                                       filtering_criteria["max_complementary_any"],
                                       rev_primer, filtering_criteria["max_diff_temp"])
    return _run_primer3(filtering_criteria["primer3_core_fpath"], arguments)


def _parse_primer3_result_for_kmer(result):
    explanation = config.EXPLANATION_LEFT_PRIMER.search(
        result).group(0).split("=")[1].rstrip()
    return explanation


def _parse_primer3_result_for_pair(result):
    explanation = config.PAIR_MSG.search(
        result).group(0).split("=")[1].rstrip()
    return explanation


def parse_primer3_results_for_kmer(results):
    explanations = []
    for result in results:
        explanations.append(
            _parse_primer3_result_for_kmer(result))
    return explanations


def parse_primer3_results_for_pair(results):
    explanations = []
    for result in results:
        explanations.append(
            _parse_primer3_result_for_pair(result))
    return explanations


def combinations_validated_by_primer3(compatibility_group_primers, primer, filtering_criteria=config.PAIR_PRIMER3_FILTERING_CRITERIA):
    for compatibility_group_primer in compatibility_group_primers:
        combination = (compatibility_group_primer, primer)
        results = _run_primer3_with_pair(combination, filtering_criteria)
        combination_explanations = parse_primer3_results_for_pair(results)
        if not all(explanation == config.OK_PRIMER3_MSG for explanation in combination_explanations):
            return False
    return True


def kmer_validated_by_primer3(kmer, filtering_criteria=config.KMER_PRIMER3_FILTERING_CRITERIA):
    results = _run_primer3_with_kmer(kmer, filtering_criteria)
    kmer_explanations = parse_primer3_results_for_kmer(results)
    if all(explanation == config.OK_PRIMER3_MSG for explanation in kmer_explanations):
        return True
    else:
        return False


def grouper(n, iterable, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return zip_longest(fillvalue=fillvalue, *args)
