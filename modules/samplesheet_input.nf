// =======================================
// MODULE: SAMPLESHEET_INPUT
// =======================================

nextflow.enable.dsl = 2

process SAMPLESHEET_INPUT {
	tag "Validating ${samplesheet.name}"
	label 'process_low'

	input:
	// Recibe path y cabecera
	tuple path(samplesheet), val(header)

	output:
	path "samplesheet_validated.csv", emit: cleaned
	path "versions.yml"         	, emit: versions // Emitir versiones

	script:
	def expected_header = 'group,replicate,input_file,fasta,gtf,technology'
	"""
	#!/bin/bash
	set -euo pipefail

	echo "Validating samplesheet header..."
	# Quitar espacios y comparar (más flexible)
	cleaned_header=\$(echo "${header}" | tr -d '[:space:]')
	cleaned_expected=\$(echo "${expected_header}" | tr -d '[:space:]')

	if [[ "\$cleaned_header" != "\$cleaned_expected" ]]; then
    	echo "ERROR: Invalid samplesheet header." >&2
    	echo "Expected: ${expected_header}" >&2
    	echo "Found:	${header}" >&2
    	exit 1
	fi
	echo "Header OK."

	# Copiar el archivo original validado
	cp "${samplesheet}" samplesheet_validated.csv

	# Verificar que el archivo copiado no está vacíoS
	if [ ! -s "samplesheet_validated.csv" ]; then
     	echo "ERROR: Failed to create or copy validated samplesheet." >&2
     	exit 1
	fi
	echo "Validated samplesheet created: samplesheet_validated.csv"

	# Versiones (bash es parte del sistema base)
	cat <<-END_VERSIONS > versions.yml
	"${task.process}":
    	bash: \$(bash --version | head -n 1)
	END_VERSIONS
	"""

	stub:
	"""
	touch samplesheet_validated.csv
	echo "${task.process}:" > versions.yml
	echo "  bash: stub" >> versions.yml
	"""
}
