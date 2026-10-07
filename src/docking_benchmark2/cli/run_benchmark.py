"""CLI for running benchmarks."""

import argparse
from pathlib import Path

from ..pipeline import run_pipeline
from ..config import load_config, load_methods_config


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Run docking benchmark pipeline 2.0",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # General arguments
    parser.add_argument('--base-dir', type=str,
                       help='Base directory for output and processed files')
    parser.add_argument('--protein-dir', type=str,
                       help='Directory containing protein PDB files')
    parser.add_argument('--ligand-dir', type=str,
                       help='Directory containing ligand CSV files')
    parser.add_argument('--methods', nargs='+',
                       choices=['qvina', 'gnina', 'plapt', 'dynamicbind', 'interformer'],
                       help='Docking methods to run')
    parser.add_argument('--protein-settings', type=str,
                       help='Path to protein preparation overrides (YAML/JSON)')
    parser.add_argument('--box-settings', type=str,
                       help='Path to box preparation overrides (YAML/JSON)')
    parser.add_argument('--interaction-config', type=str,
                       help='Path to interaction config JSON file')
    
    # Gnina arguments
    parser.add_argument('--gnina-binary', type=str,
                       help='Path to Gnina binary')
    parser.add_argument('--gnina-exhaustiveness', type=int,
                       help='Gnina exhaustiveness parameter')
    parser.add_argument('--gnina-num-modes', type=int,
                       help='Number of binding modes to generate')
    parser.add_argument('--gnina-use-cnn', dest='gnina_use_cnn',
                       action='store_true', default=None,
                       help='Use CNN scoring for Gnina (--cnn --no_gpu)')
    parser.add_argument('--gnina-no-cnn', dest='gnina_use_cnn',
                       action='store_false',
                       help='Disable CNN scoring for Gnina (--no_gpu only; AffiTox default)')
    parser.add_argument('--gnina-random-seed', type=int,
                       help='Random seed for Gnina (default: 42)')
    parser.add_argument('--gnina-docking-timeout', type=int,
                       help='Timeout per ligand in seconds for Gnina (default: 600)')
    
    # QVina arguments
    parser.add_argument('--qvina-binary', type=str,
                       help='Path to QVina binary')
    parser.add_argument('--qvina-exhaustiveness', type=int,
                       help='QVina exhaustiveness parameter')
    parser.add_argument('--qvina-num-modes', type=int,
                       help='Maximum number of binding modes to generate for QVina (default: 10)')
    parser.add_argument('--qvina-random-seed', type=int,
                       help='Random seed for QVina (for Python random/numpy, default: 42)')
    
    # PLAPT arguments
    parser.add_argument('--plapt-path', type=str,
                       help='Path to PLAPT directory')
    parser.add_argument('--plapt-device', type=str, choices=['cuda', 'cpu'],
                       help='Device to use for PLAPT (cuda or cpu)')
    parser.add_argument('--plapt-batch-size', type=int,
                       help='Batch size for PLAPT molecule encoding')
    parser.add_argument('--plapt-random-seed', type=int,
                       help='Random seed for PLAPT (for Python random/numpy/torch, default: 42)')
    
    # DynamicBind arguments
    parser.add_argument('--dynamicbind-path', type=str,
                       help='Path to DynamicBind directory')
    parser.add_argument('--dynamicbind-device', type=int,
                       help='CUDA device ID for DynamicBind')
    parser.add_argument('--dynamicbind-samples-per-complex', type=int,
                       help='Number of samples to generate per complex')
    parser.add_argument('--dynamicbind-savings-per-complex', type=int,
                       help='Number of samples to save per complex')
    parser.add_argument('--dynamicbind-inference-steps', type=int,
                       help='Number of inference steps (movie frames)')
    parser.add_argument('--dynamicbind-num-workers', type=int,
                       help='Number of workers for relaxation')
    parser.add_argument('--dynamicbind-use-relax', action='store_true',
                       help='Use relaxation for final structures')
    parser.add_argument('--dynamicbind-no-relax', dest='dynamicbind_use_relax',
                       action='store_false', help='Disable relaxation')
    parser.add_argument('--dynamicbind-rigid-protein', action='store_true',
                       help='Use rigid protein (disable protein dynamics)')
    parser.add_argument('--dynamicbind-random-seed', type=int,
                       help='Random seed for DynamicBind')
    parser.add_argument('--dynamicbind-python-env', type=str,
                       help='Path to Python executable in dynamicbind environment')
    parser.add_argument('--dynamicbind-relax-python-env', type=str,
                       help='Path to Python executable in relax environment')
    
    # Interformer arguments
    parser.add_argument('--interformer-repo-path', type=str,
                       help='Path to Interformer_bench repository')
    parser.add_argument('--interformer-conda-env', type=str,
                       help='Conda environment name for Interformer (default: interformer)')
    parser.add_argument('--interformer-gpus', type=int,
                       help='Number of GPUs for Interformer (default: 1)')
    parser.add_argument('--interformer-omp-threads', type=int,
                       help='Number of OMP threads for Interformer (default: 8)')
    parser.add_argument('--interformer-pocket-radius', type=float,
                       help='Pocket radius for Interformer in Angstroms')
    parser.add_argument('--interformer-max-ligands', type=int,
                       help='Maximum number of ligands to process per protein')
    parser.add_argument('--interformer-energy-checkpoint', type=str,
                       help='Path to Interformer energy model checkpoint')
    parser.add_argument('--interformer-affinity-checkpoint', type=str,
                       help='Path to Interformer affinity model checkpoint (glob pattern)')
    
    # Config files
    parser.add_argument('--config', type=str,
                       help='Path to config YAML file')
    parser.add_argument('--methods-config', type=str,
                       help='Path to methods config YAML file')
    
    # Stage selection - can be comma-separated list or single stage
    parser.add_argument('--stage', type=str,
                       default='all',
                       help='Which stage(s) to run: preparation, docking, aggregation, analysis, all, or comma-separated list (e.g., "docking,aggregation,analysis")')
    
    args = parser.parse_args()
    
    # Load config
    if args.config:
        config = load_config(Path(args.config))
    else:
        config = load_config()
    
    # Override with command line arguments
    if args.base_dir:
        config['base_dir'] = args.base_dir
    if args.protein_dir:
        config['protein_dir'] = args.protein_dir
    if args.ligand_dir:
        config['ligand_dir'] = args.ligand_dir
    if args.methods:
        config['methods'] = args.methods
    if args.protein_settings:
        config['protein_settings_file'] = args.protein_settings
    if args.box_settings:
        config['box_settings_file'] = args.box_settings
    if args.interaction_config:
        config['interaction_config_file'] = args.interaction_config
    
    # Load methods config
    if args.methods_config:
        methods_config = load_methods_config(Path(args.methods_config))
    else:
        methods_config = load_methods_config()
    
    # Override Gnina config with CLI args
    if args.gnina_binary:
        methods_config.setdefault('gnina', {})['binary'] = args.gnina_binary
    if args.gnina_exhaustiveness:
        methods_config.setdefault('gnina', {})['exhaustiveness'] = args.gnina_exhaustiveness
    if args.gnina_num_modes:
        methods_config.setdefault('gnina', {})['num_modes'] = args.gnina_num_modes
    if args.gnina_use_cnn is not None:
        methods_config.setdefault('gnina', {})['use_cnn'] = args.gnina_use_cnn
    if args.gnina_random_seed:
        methods_config.setdefault('gnina', {})['random_seed'] = args.gnina_random_seed
    if args.gnina_docking_timeout:
        methods_config.setdefault('gnina', {})['docking_timeout'] = args.gnina_docking_timeout
    
    # Override QVina config with CLI args
    if args.qvina_binary:
        methods_config.setdefault('qvina', {})['binary'] = args.qvina_binary
    if args.qvina_exhaustiveness:
        methods_config.setdefault('qvina', {})['exhaustiveness'] = args.qvina_exhaustiveness
    if args.qvina_num_modes:
        methods_config.setdefault('qvina', {})['num_modes'] = args.qvina_num_modes
    if args.qvina_random_seed:
        methods_config.setdefault('qvina', {})['random_seed'] = args.qvina_random_seed
    
    # Override PLAPT config with CLI args
    if args.plapt_path:
        methods_config.setdefault('plapt', {})['plapt_path'] = args.plapt_path
    if args.plapt_device:
        methods_config.setdefault('plapt', {})['device'] = args.plapt_device
    if args.plapt_batch_size:
        methods_config.setdefault('plapt', {})['batch_size'] = args.plapt_batch_size
    if args.plapt_random_seed:
        methods_config.setdefault('plapt', {})['random_seed'] = args.plapt_random_seed
    
    # Override DynamicBind config with CLI args
    if args.dynamicbind_path:
        methods_config.setdefault('dynamicbind', {})['dynamicbind_path'] = args.dynamicbind_path
    if args.dynamicbind_device is not None:
        methods_config.setdefault('dynamicbind', {})['device'] = args.dynamicbind_device
    if args.dynamicbind_samples_per_complex:
        methods_config.setdefault('dynamicbind', {})['samples_per_complex'] = args.dynamicbind_samples_per_complex
    if args.dynamicbind_savings_per_complex:
        methods_config.setdefault('dynamicbind', {})['savings_per_complex'] = args.dynamicbind_savings_per_complex
    if args.dynamicbind_inference_steps:
        methods_config.setdefault('dynamicbind', {})['inference_steps'] = args.dynamicbind_inference_steps
    if args.dynamicbind_num_workers:
        methods_config.setdefault('dynamicbind', {})['num_workers'] = args.dynamicbind_num_workers
    if hasattr(args, 'dynamicbind_use_relax'):
        methods_config.setdefault('dynamicbind', {})['use_relax'] = args.dynamicbind_use_relax
    if args.dynamicbind_rigid_protein:
        methods_config.setdefault('dynamicbind', {})['protein_dynamic'] = False
    if args.dynamicbind_random_seed:
        methods_config.setdefault('dynamicbind', {})['random_seed'] = args.dynamicbind_random_seed
    if args.dynamicbind_python_env:
        methods_config.setdefault('dynamicbind', {})['python_env'] = args.dynamicbind_python_env
    if args.dynamicbind_relax_python_env:
        methods_config.setdefault('dynamicbind', {})['relax_python_env'] = args.dynamicbind_relax_python_env
    
    # Override Interformer config with CLI args
    if args.interformer_repo_path:
        methods_config.setdefault('interformer', {})['repo_path'] = args.interformer_repo_path
    if args.interformer_conda_env:
        methods_config.setdefault('interformer', {})['conda_env'] = args.interformer_conda_env
    if args.interformer_gpus:
        methods_config.setdefault('interformer', {})['gpus'] = args.interformer_gpus
    if args.interformer_omp_threads:
        methods_config.setdefault('interformer', {})['omp_threads'] = args.interformer_omp_threads
    if args.interformer_pocket_radius:
        methods_config.setdefault('interformer', {})['pocket_radius'] = args.interformer_pocket_radius
    if args.interformer_max_ligands:
        methods_config.setdefault('interformer', {})['max_ligands'] = args.interformer_max_ligands
    if args.interformer_energy_checkpoint:
        methods_config.setdefault('interformer', {})['energy_checkpoint'] = args.interformer_energy_checkpoint
    if args.interformer_affinity_checkpoint:
        methods_config.setdefault('interformer', {})['affinity_checkpoint'] = args.interformer_affinity_checkpoint
    
    # Run pipeline
    run_pipeline(
        config=config,
        methods_config=methods_config,
        stage=args.stage,
    )


if __name__ == '__main__':
    main()

