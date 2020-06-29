all:
	pandoc --template projeto.latex --filter pandoc-fignos --from markdown+multiline_tables --output simulation_methods.pdf simulation_methods.md 