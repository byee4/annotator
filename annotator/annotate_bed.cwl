#!/usr/bin/env python

cwlVersion: v1.0

class: CommandLineTool

baseCommand: [annotate-bed]

inputs:

  input:
    type: File
    inputBinding:
      position: 1
      prefix: --input
    label: "input BED6 file"
    doc: "input unannotated BED6 file"

  output:
    type: string
    inputBinding:
      position: 2
      prefix: --output
    label: "output tsv file"
    doc: "annotated tabbed file"

  gtfdb:
    type: File
    inputBinding:
      position: 3
      prefix: --gtfdb

  transcriptPriorityFile:
    type: File
    inputBinding:
      position: 4
      prefix: --transcript-priority-file

  genePriorityFile:
    type: File
    inputBinding:
      position: 5
      prefix: --gene-priority-file

  unstranded:
    type: boolean
    default: false
    inputBinding:
      position: 6
      prefix: --unstranded

  appendChr:
    type: boolean
    default: False
    inputBinding:
      position: 7
      prefix: --append-chr

  species:
    type: string
    inputBinding:
      position: 8
      prefix: --species

  # limitChromsTo:
  #   type: string[]
  #   inputBinding:
  #     position: 9

outputs:

  output_file:
    type: File
    outputBinding:
      glob: $(inputs.output)
    label: "output"
    doc: "File containing output of the annotation program"


