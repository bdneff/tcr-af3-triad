// hla2_rmsd.nf -- CDR RMSD for the five class II triads folded by hla2.nf.
//
// Copy to $WS/tcrtrifold-experiments/workflows/hla2_rmsd.nf and run from there.
//
//   nextflow run ./workflows/hla2_rmsd.nf \
//       --input_triad $PWD/data/hla2/triad/staged/hla2_triad.segid.parquet \
//       --triad_inf_dir $PWD/data/hla2/triad/inference \
//       -profile gemini -resume
//
// Same four stages pdb.nf uses for its own RMSDs, in the same order, with the same
// processes -- so `cdr_rmsd_af3_0` here means exactly what it means for the published
// points on the panel:
//
//   FORMAT_TRUE_PDBS          cut the deposited crystal into canonical A/B/C/D/E
//   TCRDOCK_GEOM_FROM_PDB     TCRdock frame for the crystal
//   TCRDOCK_GEOM_FROM_AF3_INFERENCE   TCRdock frame for the prediction
//   COMPUTE_RMSD_AF3          superpose on the MHC frame, IMGT-number, CA RMSD over CDRs
//
// Input must be the *segid* parquet: FORMAT_TRUE_PDBS selects `segid <X>` from the
// deposit, and the staged table carries role labels ("alpha"/"beta") rather than chain
// IDs until 02_derive_segids.py fills them in.
//
// No GPU -- this reads predictions that already exist.
//
// RETRIEVING THE RESULT. Do not trust a green summary; hla2.nf reported every process at
// 100% and "Succeeded: 67" while writing nothing to the path it was given. Find the
// output parquet through the log instead:
//
//   grep "COMPUTE_RMSD_AF3" logs/hla2_rmsd/.nextflow.log \
//     | grep -oE "workDir: [^]]+" | sed 's/workDir: //' | sort -u
//   # then copy *.parquet out of that directory
//
// and confirm it has five rows before using it.

include { FORMAT_TRUE_PDBS } from './subworkflows/local/cleaning'
include { TCRDOCK_GEOM_FROM_PDB;
          TCRDOCK_GEOM_FROM_INFERENCE as TCRDOCK_GEOM_FROM_AF3_INFERENCE;
          COMPUTE_RMSD as COMPUTE_RMSD_AF3 } from './subworkflows/local/extract_feat'

workflow {

    main:

    triad = Channel.fromPath(params.input_triad)
    triad_inf_dir = file(params.triad_inf_dir).toUriString()

    cleaned_pdbfiles = FORMAT_TRUE_PDBS(triad)

    triad_tcrdock_pdb = TCRDOCK_GEOM_FROM_PDB(triad, cleaned_pdbfiles)
    triad_tcrdock_af3 = TCRDOCK_GEOM_FROM_AF3_INFERENCE(triad, triad_inf_dir,
                                                        Channel.value("af3"))

    COMPUTE_RMSD_AF3(triad, triad_tcrdock_pdb, triad_tcrdock_af3, cleaned_pdbfiles,
                     triad_inf_dir, Channel.value("af3"))

    publish:
    triad_rmsd_af3 = COMPUTE_RMSD_AF3.out.triad
    cleaned_pdbfiles = cleaned_pdbfiles
}
