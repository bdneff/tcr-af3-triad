// hla2.nf -- fold the class II triads the published pdb.nf run never covered.
//
// Copy to $WS/tcrtrifold-experiments/workflows/hla2.nf and run from there; the
// include paths are relative to that directory, and the repo keeps one .nf per input
// set (pdb, cresta, iedb_I, iedb_II, cross_reactivity, ...) so this follows the
// convention rather than patching pdb.nf.
//
// Why a separate entry point: pdb.nf's workflow is unnamed (so nothing can include it)
// and starts from DISCOVERY inputs -- the Table S1 benchmark CSV and STCRDab's
// db_summary.dat -- which it feeds to CLEAN_PDB to resolve triads from RCSB. Our five
// triads are already resolved, and appear in neither of those inputs. So we enter at
// the triad-parquet stage instead, which is exactly what the af3_adapter subworkflows
// take.
//
// This deliberately stops after inference. The metrics come from separate steps:
// COMPUTE_RMSD needs cleaned crystals from FORMAT_TRUE_PDBS (which needs the *_segid
// columns our staged table does not yet carry), and PTI-PAE comes from the standalone
// workflows/bin/extract_triad_conf_feat.py. Inference is the only part that needs a
// GPU and the only part worth queueing now.
//
//   nextflow run ./workflows/hla2.nf \
//       --input_triad data/hla2/triad/staged/hla2_triad.staged.parquet \
//       --triad_inf_dir $PWD/data/hla2/triad/inference \
//       -profile gemini -resume
//
// Pass --triad_inf_dir as an ABSOLUTE path, and check the count when it finishes:
//     ls data/hla2/triad/inference | wc -l
// A green summary is not evidence the predictions landed -- see the note below.
//
// MSAs are ON and there is no --norun_data_pipeline anywhere, matching the published
// predictions: the manuscript Methods specify MSAs on, templates on, seeds 1-5, one
// diffusion sample, 10 recycles, and MSA_WORKFLOW:RUN_MSA appears in the published run
// log. Folding these from bare sequence would not be comparable to the other points on
// the panel.

include { MSA_FROM_TRIAD_PARQUET;
          UNBATCHED_INFERENCE_FROM_TRIAD_PARQUET as UNBATCHED_INFERENCE_HLA2;
          NOOP_DEP as DEPEND_HLA2_ON_INFERENCE } from './subworkflows/local/af3_adapter'

workflow {

    main:

    triad = Channel.fromPath(params.input_triad)

    // Take the inference directory EXPLICITLY. Deriving it from $workflow.outputDir
    // silently failed: -output-dir only drives the `output {}` publishing mechanism
    // (the preview feature Nextflow warns about), which this file does not define, so
    // the path never resolved. Every process reported 100% and Succeeded: 67 while
    // data/hla2/triad/inference did not exist -- AF3 runs with `--output_dir=.` and had
    // written into its task directories under /scratch/$USER/work, which is swept.
    // Nothing was lost, but only counting the outputs revealed it.
    triad_inf_dir = file(params.triad_inf_dir).toUriString()

    MSA_FROM_TRIAD_PARQUET(triad)
    triad_msa_done_token = MSA_FROM_TRIAD_PARQUET.out.new_msa_list

    UNBATCHED_INFERENCE_HLA2(triad, triad_inf_dir, triad_msa_done_token)
    triad_meta_inf = UNBATCHED_INFERENCE_HLA2.out.new_meta_inf

    triad_folded = DEPEND_HLA2_ON_INFERENCE(triad, triad_meta_inf.toList())

    publish:
    triad_cleaned = triad_folded
    triad_meta_inf = triad_meta_inf
}
