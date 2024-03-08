all:
	#pandoc --template projeto.latex --filter pandoc-xnos --from markdown+multiline_tables --pdf-engine=xelatex --output simulation_methods.pdf simulation_methods.md 
	pandoc --filter pandoc-crossref --citeproc --pdf-engine=xelatex --output simulation_methods.pdf simulation_methods.md 
	
	#pandoc --template projeto.latex --filter pandoc-xnos --from markdown+multiline_tables --output simulation_methods.docx simulation_methods.md 
	#cp simulation_methods.docx ~/gdrive