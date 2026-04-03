#!/bin/bash

BASEDIR="/ceph/cms/store/user/tvami/DiMuonPlusX/Ntuples_v1.0.1/Data"
UNMERGEDBASE="/ceph/cms/store/user/tvami/DiMuonPlusX/Ntuples_v1.0.1/unmerged"
# BASEDIR="/ceph/cms/store/user/tvami/DiMuonPlusX/Ntuples_v1.0.1_trigger_study/Data"
# UNMERGEDBASE="/ceph/cms/store/user/tvami/DiMuonPlusX/Ntuples_v1.0.1_trigger_study/unmerged"
BATCH_SIZE=50

for region in sr vr; do
    for era in E F G H I; do
        inputdir="${BASEDIR}/${region}"
        output="${inputdir}/skimmed_${region}_ScoutingPFRun3_Run2024${era}.root"
        mapfile -t files < <(ls ${inputdir}/skimmed_${region}_ScoutingPFRun3_Run2024${era}_*_job*.root 2>/dev/null)
        # output="${inputdir}/trigger_study_${region}_ScoutingPFRun3_Run2024${era}.root"
        # mapfile -t files < <(ls ${inputdir}/trigger_study_${region}_ScoutingPFRun3_Run2024${era}_*_job*.root 2>/dev/null)
        tmp_existing=""

        # Filter out corrupt ROOT files
        good_files=()
        for f in "${files[@]}"; do
            if python3 -c "import ROOT; f=ROOT.TFile.Open('${f}'); exit(0 if f and not f.IsZombie() and f.GetListOfKeys() else 1)" 2>/dev/null; then
                good_files+=("$f")
            else
                echo "WARNING: Skipping corrupt file: $f"
            fi
        done
        files=("${good_files[@]}")

        nfiles=${#files[@]}
        if [ "$nfiles" -eq 0 ]; then
            echo "No files found for ${region} Run2024${era}, skipping."
            continue
        fi

        UNMERGED="${UNMERGEDBASE}/Data_${region}"
        mkdir -p "${UNMERGED}"

        if [ -f "${output}" ]; then
            echo "Existing merged file found, will merge ${nfiles} new files into it."
            tmp_existing="${output%.root}_old.root"
            mv "${output}" "${tmp_existing}"
            files=("${tmp_existing}" "${files[@]}")
            nfiles=${#files[@]}
        fi

        echo "Merging ${nfiles} files for ${region} Run2024${era} -> ${output}"

        if [ "$nfiles" -le "$BATCH_SIZE" ]; then
            hadd -f "${output}" "${files[@]}"
        else
            tmpdir=$(mktemp -d "${inputdir}/hadd_tmp_${era}_XXXX")
            batch=0
            failed=0
            for ((i=0; i<nfiles; i+=BATCH_SIZE)); do
                batch_files=("${files[@]:i:BATCH_SIZE}")
                batch_out="${tmpdir}/batch_${batch}.root"
                echo "  Batch ${batch}: merging ${#batch_files[@]} files"
                hadd -f "${batch_out}" "${batch_files[@]}"
                if [ $? -ne 0 ]; then
                    echo "ERROR: hadd failed on batch ${batch} for ${region} Run2024${era}"
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
                    echo "ERROR: final hadd merge failed for ${region} Run2024${era}"
                fi
            fi

            rm -rf "${tmpdir}"
        fi

        if [ $? -eq 0 ] && [ -f "${output}" ]; then
            echo "Successfully created ${output}"
            [ -n "${tmp_existing}" ] && [ -f "${tmp_existing}" ] && rm "${tmp_existing}"
            mapfile -t job_files < <(printf '%s\n' "${files[@]}" | grep '_job')
            if [ ${#job_files[@]} -gt 0 ]; then
                echo "Moving ${#job_files[@]} input job files to ${UNMERGED}/"
                mv "${job_files[@]}" "${UNMERGED}/"
            fi
        else
            echo "ERROR: hadd failed for ${region} Run2024${era}, keeping original files in place."
            [ -n "${tmp_existing}" ] && [ -f "${tmp_existing}" ] && mv "${tmp_existing}" "${output}"
        fi
        echo ""
    done
done
