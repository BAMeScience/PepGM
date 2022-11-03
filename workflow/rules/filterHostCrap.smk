rule AddHostandCrap:
    input:
       DatabaseDirectory+'crap.fasta',
       DatabaseDirectory+HostName +'.fasta'

    output: DatabaseDirectory+ HostName+'+crap.fasta'

    shell:'cat {input} > {output}'  

rule RemoveDuplicatesHostandCrap:
    input: DatabaseDirectory+ HostName+'+crap.fasta'
    output: DatabaseDirectory+ HostName+'+crap_UNI.fasta'
    conda: 'envs/graphenv.yml'
    shell:  'seqkit rmdup -s {input} > {output}'

rule AddDecoysHostandCrap:
    input: DatabaseDirectory+HostName+'+crap_UNI.fasta'
    output: 
          DatabaseDirectory+HostName+'+crap_UNI_concatenated_target_decoy.fasta',
          DatabaseDirectory+HostName+'+crap_UNI_decoy.fasta'
    conda: 'envs/decoy.yml'
    shell: 'decoypyrat {input} -o {output[1]} -d DECOY && cat {input} {output[1]} > {output[0]}'


### X!Tandem ###

XTandemOutputFilter = XTandemDir + 'tandem_HostCrap.output.xml'


# create input file for X!Tandem
rule CreateXTandemInputFilter:
    input:
        spectra = DataDirectory+SampleName+'/'+SampleName+SpectraFileType,
        taxonomy = XTandemDir + 'tandem_Filter.taxonomy.xml'
    output:
        XTandemDir+'tandem_Filter.input.xml'
    params:
        xtandem_output = XTandemOutputFilter
    run:
        f_out = open(output[0], "w")

        lines = ["<?xml version=\"1.0\"?>",
                "<bioml>"]

        lines.append(header_note("list path parameters"))
        lines.append(input_note("list path, default parameters", XTANDEM_DEFAULT))
        lines.append(input_note("list path, taxonomy information", input.taxonomy))

        lines.append(header_note("spectrum parameters"))
        lines.append(input_note("spectrum, path", input.spectra))

        lines.append(input_note("spectrum, fragment monoisotopic mass error", XTANDEM_FMME))
        lines.append(input_note("spectrum, fragment monoisotopic mass error units", XTANDEM_FMMEU))

        lines.append(input_note("spectrum, parent monoisotopic mass error plus", XTANDEM_PMMEP))
        lines.append(input_note("spectrum, parent monoisotopic mass error minus", XTANDEM_PMMEM))
        lines.append(input_note("spectrum, parent monoisotopic mass error units", XTANDEM_PMMEU))

        lines.append(header_note("spectrum conditioning parameter"))
        lines.append(input_note("spectrum, threads", threads))

        if XTANDEM_MODS_FIX or XTANDEM_MODS_VAR:
            lines.append(header_note("residue modification parameters"))
        if XTANDEM_MODS_FIX:
            lines.append(input_note("residue, modification mass", XTANDEM_MODS_FIX))
        if XTANDEM_MODS_VAR:
            lines.append(input_note("residue, potential modification mass", XTANDEM_MODS_VAR) )
        if XTANDEM_MODS_VAR_NTERM:
            lines.append(header_note("model refinement parameters"))
            lines.append(input_note("refine", "yes"))
            lines.append(input_note("refine, potential N-terminus modifications", XTANDEM_MODS_VAR_NTERM))

        lines.append(header_note("protein parameters"))
        lines.append(input_note("protein, taxon", "all_viruses"))
        # enzyme is trypsin by default

        lines.append(header_note("scoring parameters"))
        lines.append(input_note("scoring, include reverse", "no"))

        lines.append(header_note("output parameters"))
        lines.append(input_note("output, path", params.xtandem_output))
        lines.append(input_note("output, path hashing", "no"))  # no date/time tag in output file name
        lines.append(input_note("output, message", ""))

        if XTANDEM_PARAS:
            lines.append(header_note("additional parameters"))
            lines.extend([input_note(param, value) for param, value in XTANDEM_PARAS.items()])

        lines.append("</bioml>")

        f_out.writelines([line + "\n" for line in lines])
        f_out.close()

#create XTandem Taxonomy file
rule XTandemTaxonomyInputFilter:
    input:
        DatabaseDirectory+HostName+'+crap_UNI_concatenated_target_decoy.fasta'
    output:
        XTandemDir+'tandem_Filter.taxonomy.xml'
    run:
        lines = ["<?xml version=\"1.0\"?>",
        "<bioml label=\"x! taxon-to-file matching list\">",
        "<taxon label=\"all_viruses\">",
        "<file format=\"peptide\" URL=\"" + DatabaseDirectory+HostName+'+crap_UNI_concatenated_target_decoy.fasta' + "\" />"]
        lines.append("</taxon>")
        lines.append("</bioml>")
        f_out = open(output[0], "w")
        f_out.writelines([line + "\n" for line in lines])
        f_out.close()

# locate default X!Tandem executable
# look for xtandem (conda) first, then tandem.exe (manual installation)
def xtandem_executable():
    import shutil
    if shutil.which("xtandem"):
        return("xtandem")
    elif shutil.which("tandem.exe"):
        return ("tandem.exe")
    else:
        print("neither xtandem or tandem.exe is available", file=sys.stderr)
        sys.exit(1)


rule ExecuteXTandemFilter:
    input:
        spectra = DataDirectory+SampleName+'/'+SampleName+SpectraFileType,
        paras = XTandemDir+'tandem_Filter.input.xml',
        taxon = XTandemDir+'tandem_Filter.taxonomy.xml'
    output:
        XTandemOutputFilter
    conda: 'envs/graphenv.yml'
    #threads:
     #   THREADS_MAX
    #benchmark:
    #    "logs/xtandem.txt"
    log:
        XTandemDir + 'xtandem.log'
    #params:
    #   exec = xtandem_executable()
    shell:
        'xtandem {input.paras} 2>&1 | tee {log}'

