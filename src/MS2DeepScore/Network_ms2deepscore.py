import argparse
from matchms import calculate_scores
from matchms.importing import load_from_mgf
from ms2deepscore import MS2DeepScore
from ms2deepscore.models import load_model
import csv
from matchms.filtering import select_by_mz, default_filters, add_fingerprint
from matchms.filtering import normalize_intensities, add_parent_mass, derive_adduct_from_name
from matchms.filtering import require_minimum_number_of_peaks

from matchms.filtering.metadata_processing.repair_not_matching_annotation import repair_not_matching_annotation
from matchms.filtering.metadata_processing.derive_inchi_from_smiles import derive_inchi_from_smiles
from matchms.filtering.metadata_processing.derive_inchi_from_smiles  import derive_inchi_from_smiles
from matchms.filtering.metadata_processing.derive_inchikey_from_inchi import derive_inchikey_from_inchi
from matchms.filtering import reduce_to_number_of_peaks
from matchms.filtering import select_by_relative_intensity
from matchms import Spectrum


def process_peaks(s):
    s = select_by_relative_intensity(s, intensity_from=0.001)
    s = reduce_to_number_of_peaks(s, n_max=1000)
    return s

def minimal_processing(spectrum):
    spectrum = default_filters(spectrum)
    spectrum = derive_adduct_from_name(spectrum)
    # spectrum = repair_not_matching_annotation(spectrum)
    spectrum = add_parent_mass(spectrum)
    spectrum = add_fingerprint(spectrum)
    # If a fingerprint couldn't be generated, discard spectrum
    if spectrum is None:
        return  None
    spectrum = normalize_intensities(spectrum)
    spectrum = select_by_mz(spectrum, mz_from=10.0, mz_to=1000.0)
    spectrum = select_by_relative_intensity(spectrum, intensity_from=0.001)
    spectrum = require_minimum_number_of_peaks(spectrum, n_required=1)
    return spectrum
def convert_extension(input_string,threshold, new_extension):
    # Split the input string into name and extension
    name, old_extension = input_string.rsplit('.', 1)
    name = name + "_" + str(threshold)

    # Check if the old extension is "mgf"
    if old_extension.lower() == 'mgf':
        # Concatenate the name with the new extension
        output_string = f"{name}.{new_extension}"
        return output_string
    else:
        # If the input string doesn't end with ".mgf", return an error message or handle it as needed
        return "Invalid input file format"
def main():
    # pass arguments
    parser = argparse.ArgumentParser(description='Using ms2deepscore model to construct the network barebone')
    parser.add_argument('-m', type=str, required=True, default="spec.mgf", help='input mgf filename')
    parser.add_argument('-w', type=str, required=True, default="MS2DeepScore_allGNPSpositive_10k_500_500_200.hdf5", help='model weights')
    args = parser.parse_args()
    input_mgf = args.m
    model = args.w
    spectra = list(load_from_mgf(input_mgf))
    spectra = [minimal_processing(x) for x in spectra]
    spectra = [process_peaks(x) for x in spectra]
    spectra = [x for x in spectra if x is not None]
    for spectrum in spectra:
        if spectrum is None:
            print("None of the spectra")
    model = load_model(model)
    similarity_measure = MS2DeepScore(model)
    scores = calculate_scores(spectra, spectra, similarity_measure, is_symmetric=True)
    spec_scores = scores.scores.data.reshape(len(spectra),len(spectra))
    threshold_list = [0.5,0.6,0.7,0.8,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,0.991,0.992,0.993,0.994,0.995,0.996,0.997,0.998,0.999,0.9995,0.9996,0.9997,0.9998,0.9999,0.99999]
    for threshold in threshold_list:
        output_file = convert_extension(input_mgf, threshold, "tsv")
        with open("../../results/MS2DeepScoreNetwork/"+output_file, 'w', newline='') as csvfile:
            # Create a CSV writer with tab as the delimiter
            csv_writer = csv.writer(csvfile, delimiter='\t')

            # Write the header
            csv_writer.writerow(['CLUSTERID1', 'CLUSTERID2', 'Cosine'])

            # Iterate through the pairwise_scores and write rows where the score is over 0.7
            for cluster_id1, row in enumerate(spec_scores):
                for cluster_id2, score in enumerate(row):
                    if cluster_id1 < cluster_id2 and score[0] > threshold:
                        csv_writer.writerow([cluster_id1+1, cluster_id2+1, score[0]])


if __name__ == '__main__':
    main()