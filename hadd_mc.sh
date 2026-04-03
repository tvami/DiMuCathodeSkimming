#!/bin/bash

#BASEDIR="/ceph/cms/store/user/tvami/DiMuonPlusX/Ntuples_v1.0.1"
BASEDIR="/ceph/cms/store/user/tvami/DiMuonPlusX/Ntuples_v1.0.1_trigger_study"
UNMERGEDBASE="${BASEDIR}/unmerged"
PREFIX="trigger_study"  # was: skimmed

BATCH_SIZE=50

# for mctype in BkgMC Signal; do
for mctype in BkgMC; do
    for region in sr vr; do  # was: sr vr
        inputdir="${BASEDIR}/${mctype}/${region}"

        # Find all unique sample names by extracting pattern: ${PREFIX}_${region}_<sample>_job*.root
        # Get unique sample names (everything between ${PREFIX}_${region}_ and _job)
        mapfile -t samples < <(ls "${inputdir}"/${PREFIX}_${region}_*_job*.root 2>/dev/null | \
            sed -E "s|.*/${PREFIX}_${region}_(.*)_job[0-9]+\.root|\1|" | sort -u)

        if [ ${#samples[@]} -eq 0 ]; then
            echo "No MC samples found in ${inputdir}, skipping."
            continue
        fi

        for sample in "${samples[@]}"; do
            output="${inputdir}/${PREFIX}_${region}_${sample}.root"
            mapfile -t files < <(ls "${inputdir}"/${PREFIX}_${region}_${sample}_job*.root 2>/dev/null)

            nfiles=${#files[@]}
            if [ "$nfiles" -eq 0 ]; then
                echo "No files found for ${region} ${sample}, skipping."
                continue
            fi

            UNMERGED="${UNMERGEDBASE}/${mctype}_${region}"
            mkdir -p "${UNMERGED}"

            echo "Merging ${nfiles} files for ${mctype} ${region} ${sample} -> ${output}"

            if [ "$nfiles" -le "$BATCH_SIZE" ]; then
                hadd -f "${output}" "${files[@]}"
            else
                tmpdir=$(mktemp -d "${inputdir}/hadd_tmp_${sample}_XXXX")
                batch=0
                failed=0
                for ((i=0; i<nfiles; i+=BATCH_SIZE)); do
                    batch_files=("${files[@]:i:BATCH_SIZE}")
                    batch_out="${tmpdir}/batch_${batch}.root"
                    echo "  Batch ${batch}: merging ${#batch_files[@]} files"
                    hadd -f "${batch_out}" "${batch_files[@]}"
                    if [ $? -ne 0 ]; then
                        echo "ERROR: hadd failed on batch ${batch} for ${region} ${sample}"
                        failed=1
                        break
                    fi
                    ((batch++))
                done

                if [ "$failed" -eq 0 ]; then
                    echo "  Final merge of ${batch} batches"
                    hadd -f "${output}" "${tmpdir}"/batch_*.root
                    if [ $? -ne 0 ]; then
                        failed=1
                        echo "ERROR: final hadd merge failed for ${region} ${sample}"
                    fi
                fi

                rm -rf "${tmpdir}"
            fi

            if [ $? -eq 0 ] && [ -f "${output}" ]; then
                echo "Successfully created ${output}"
                echo "Moving ${nfiles} input files to ${UNMERGED}/"
                mv "${files[@]}" "${UNMERGED}/"
            else
                echo "ERROR: hadd failed for ${region} ${sample}, keeping original files in place."
            fi
            echo ""
        done
    done
done
