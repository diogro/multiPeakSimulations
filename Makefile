all:
	pandoc --template projeto.latex --filter pandoc-fignos --output simulation_methods.pdf simulation_methods.md 