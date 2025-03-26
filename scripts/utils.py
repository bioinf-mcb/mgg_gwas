import pandas as pd
from Bio import SeqIO
from pathlib import Path
import os
import sys

def create_symlink(src_abs_path, dst_abs_path):
    """
    Create a relative symbolic link from absolute paths.
    
    Parameters:
    src_abs_path (str): The absolute path to the source file or directory.
    dst_abs_path (str): The absolute path to the destination where the symlink will be created.
    """

    # Convert the input paths to Path objects and resolve them to absolute paths
    src_path = Path(src_abs_path)
    dst_path = Path(dst_abs_path)

    # Get the directory where the symlink will be created
    symlink_dir = dst_path.parent

    # Calculate the relative path from the symlink's directory to the source
    relative_src_path = os.path.relpath(src_path, symlink_dir)
    
    # is symlink
    if dst_path.is_symlink():
        dst_path.unlink()

    # create symlink
    try:
        os.symlink(relative_src_path, dst_path)
    except OSError as e:
        print(f'Error creating symbolic link: {e}')


def symlink_directories(src_dir, dst_dir, file_pattern='*', src_paths=None, retain_structure=False):
    """
    Process files from a source directory and map them to a destination directory
    by creating relative symbolic links. Optionally, retain the folder structure.

    Parameters:
    - src_dir: Path to the source directory where files are located.
    - dst_dir: Path to the destination directory where symlinks will be created.
    - file_pattern: Pattern to match files in the source directory (e.g., '*.fasta').
    - src_paths: Optional list of source paths to use. If not provided, the function
                 will find files based on the file_pattern.
    - retain_structure: If True, retains the folder structure of src_dir in dst_dir.

    Returns:
    - Tuple of (list of source file paths, list of destination file paths as strings).
    """
    src_dir = Path(src_dir).resolve()
    dst_dir = Path(dst_dir).resolve()

    # If src_paths is not provided, gather files using a recursive glob
    if src_paths is None:
        src_paths = list(map(str, src_dir.rglob(file_pattern)))
    else:
        src_paths = list(map(str, src_paths))

    dst_paths = []
    for src in src_paths:
        src_path = Path(src)
        if retain_structure:
            # Compute the path relative to the source directory
            try:
                relative_path = src_path.relative_to(src_dir)
            except ValueError:
                raise ValueError(f"Source file {src_path} is not under the source directory {src_dir}")
            dst_path = dst_dir / relative_path
        else:
            dst_path = dst_dir / src_path.name
        dst_paths.append(str(dst_path))

    # Ensure the destination directory and all necessary parent directories exist
    Path(dst_dir).mkdir(exist_ok=True, parents=True)
    for dst in dst_paths:
        Path(dst).parent.mkdir(exist_ok=True, parents=True)

    # Create symlinks
    for src, dst in zip(src_paths, dst_paths):
        create_symlink(src, dst)

    return (src_paths, dst_paths)


def params_and_folder_compability(base_dir, params):
    """ check if the folder name entered manually is compatible with params """

    # params in folder name
    folder_name_manual = Path(base_dir).stem
    folder_parents_manual = Path(base_dir).parents
    folder_name_manual = '_'.join(folder_name_manual.split('_')[-3:])

    # actual params
    nbootstrap_actual = int(params['gwas']['nbootstrap'])
    pc_pindentity_actual, pc_coverage_actual = (str(int(float(params['mmseqs']['min_identity']) * 100))).zfill(2), (str(int(float(params['mmseqs']['min_coverage']) * 100))).zfill(2)
    annotation_tool_actual = str(params['input']['annotation_tool']).upper()
    
    folder_name_actual = f'BTSP{nbootstrap_actual}x_PCI{pc_pindentity_actual}C{pc_coverage_actual}_{annotation_tool_actual}'
    
    if folder_name_manual == folder_name_actual: pass
    else: 
        print('ERROR! Folder name is incorrect: ')
        print(f'Incorrect folder name: {folder_name_manual}')
        print(f'Correct folder name: {folder_name_actual}')

        sys.exit("Program terminated")


def get_alignments_manual(proteins_fasta, mmseqs_clusters, pcs, alignments_dir='/Users/januszkoszucki/Downloads/alignments_dir'):

    # create
    Path(alignments_dir).mkdir(exist_ok=True, parents=True)

    # load
    clusters_df = pd.read_csv(mmseqs_clusters, sep='\t')

    # filter
    filt_pcs = clusters_df['PC'].isin(pcs)
    clusters_df = clusters_df.loc[filt_pcs].reset_index(drop=True)

    # each PC
    for pc, group in clusters_df.groupby('PC'):

        # path
        fasta_file = Path(alignments_dir, f'{pc}.fasta')

        # checkpoint
        if Path(fasta_file).exists(): continue

        # alignment
        fasta = []
        for row in group.itertuples():
            for r in SeqIO.parse(proteins_fasta, 'fasta'):
                if row.proteinID == r.id:
                    fasta.append(f'>{row.proteinID}\n{r.seq}\n')
                    break
        fasta = ''.join(fasta)

        # save
        with open(fasta_file, 'w') as f:
            f.write(fasta)
