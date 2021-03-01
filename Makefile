# several sample targets

NOTES_src  = NOTES.md
NOTES_base = $(basename $(NOTES_src))

debug:
	mkdir -p build_debug && cd build_debug && cmake .. -DCMAKE_BUILD_TYPE=Debug -DCUDA_ARCH="35 52 60" -G Ninja && ninja -j4

tags:
	ctags -R --c++-kinds=+p --fields=+iaS --extra=+q .

notes:
	pandoc --pdf-engine=xelatex -V mainfont="PT Sans" $(NOTES_src) -o $(NOTES_base).pdf
	pandoc --include-in-header=fontoptions_with_slide_number.tex -s -t beamer -V theme:Warsaw --highlight-style pygments  --pdf-engine=xelatex -V fontsize=10pt $(NOTES_src) -o $(NOTES_base)_preso.pdf
