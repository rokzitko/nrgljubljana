SITE ?= auger:/var/www/html/nrgljubljana/

all: site

site: mkdocs.yml README.mkdocs serve $(shell find docs -type f)
	python3 -m mkdocs build -f mkdocs.yml

upload: site
	scp -r site/* ${SITE}

serve:
	python3 -m mkdocs serve -f mkdocs.yml

clean:
	rm -rf site/
