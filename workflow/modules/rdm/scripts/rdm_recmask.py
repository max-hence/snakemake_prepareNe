# Just randomly subsample the rec map as a control
from pandas import read_csv
import numpy as np
from argparse import ArgumentParser

def rdm_recmask(map_path:str, output_bed:str):
    """Keep randmoly a third of the map

    Args:
        recmap (str): _description_
        output_bed (str): _description_
    """
    
    recmap_df = read_csv(map_path, delimiter=" ")
    recmap_df["start"] = recmap_df["start"].astype('int')
    recmap_df["end"] = recmap_df["end"].astype('int')
    mask = np.zeros(len(recmap_df), dtype=bool)
    mask[:len(recmap_df) // 3] = True

    # Mélanger les indices des lignes où True
    np.random.shuffle(mask)

    # Appliquer le masque
    sample_df = recmap_df[mask]

    # Sauvegarder le résultat dans un nouveau fichier
    sample_df.to_csv(output_bed, sep=' ', index=False)


if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('-i', '--input', type=str,
        help="Specifies the parameters list"
    )
    parser.add_argument('-o', '--output', type=str,
        help="Specifies the output path"
    )
    args = parser.parse_args()

    rdm_recmask(args.input, args.output)