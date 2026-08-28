SITE ?= auger:/var/www/html/nrgljubljana/site/

.PHONY: all upload serve clean

all: site

site: mkdocs.yml README.mkdocs $(shell find docs -type f)
	python3 -m mkdocs build --strict -f mkdocs.yml

upload: site
	scp -O -r site/* ${SITE}

serve:
	python3 -m mkdocs serve -f mkdocs.yml

clean:
	rm -rf site/
