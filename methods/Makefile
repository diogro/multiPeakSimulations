all:
	#pandoc --template projeto.latex --filter pandoc-xnos --from markdown+multiline_tables --pdf-engine=xelatex --output simulation_methods.pdf simulation_methods.md 
	pandoc --filter pandoc-crossref --citeproc --pdf-engine=xelatex --output simulation_methods.pdf simulation_methods.md 
	
docx:
	pandoc --filter pandoc-crossref --citeproc --output simulation_methods.docx simulation_methods.md 
	