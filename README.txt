Stony Brook University PhD thesis template for LaTeX
----------------------------------------------------

Benjamin Hornberger, May 2007

This is the set of files I used to prepare my thesis in May 2007. It
is not an official (university-provided) template, so it comes with
no warranty whatsoever. It just worked for me, and grad school
accepted it. The origin of the class file is unclear, I just made
some changes to make the front matter look a little bit nicer (to my
taste) and to fix some bugs with the table of contents. I know that
there are more templates around; feel free to use what you like
best.

Make sure to get the grad schools' thesis guidelines
(www.grad.stonybrook.edu -> Graduation Information -> Guidelines to
the preparation of theses and dissertations) and verify yourself
that everything is fine in the end.

I used MiKTeX 2.5 on Windows XP and the WinEdt 5.4 editor to prepare
my thesis, but in principle this set of files should work with any
LaTeX system on any platform. I used latex and dvipdfm to produce
the final PDF file, but it should work with dvips/ps2pdf, or
pdflatex as well (in which case you should change the driver in the
hyperref.cfg file).

Beyond these few instructions, you will obviously need some general
knowledge about LaTeX as such, and your TeX system in particular. I
can recommend:
* Kopka / Daly: Guide to LaTeX (Book)
* The not so short introduction to LaTeX 2e (on the web; search
  Google)
* The comp.text.tex newsgroup (search the archives at
  groups.google.com for specific questions -- most stuff has been
  asked and answered in the past)
* Documentation to the packages used (should come with your LaTeX
  installation, or on the web)
* MiKTeX documentation and dojo (miktex.org, dojo.miktex.org) for
  MiKTeX-specific issues, default configuration files, paths to
  fonts and BiBTeX files, etc.)
* winedt.com, winedt.org for WinEdt-specific issues (how to set up
  and make best use of projects, how WinEdt knows about your BiBTeX
  files etc.)

This template package consists of the following files:

* usbthesis.cls
General document class file. Keep in the same directory as
everything else, or in the TeX search path. You shouldn't have to
make any changes to it, unless you want to change the appearance of
the front matter or you find further bugs beyond the ones I found.
To produce the extra abstract printout with the advisor's name on it
(see the grad school's guidelines), you can uncomment line 179.

* thesis.tex
This is your thesis "root" file. It loads all the packages which I
needed / found useful -- edit to your heart's content. In the
bibliography section, there are some "%input" directives which help
WinEdt find your BiBTeX files for automatic completion of \cite{}
commands. This only works with the WinEdt editor. To the TeX system,
this is just a comment.

* mathdefs.tex
This file defines a number of special commands for easier and more
consistent math typesetting (like units etc.). Most are taken from
the article which is mentioned at the top of the file. Note that
also the \vec{} command is redefined to use a boldface letter rather
than an arrow; if you don't like it, comments out that part. I also
recommend using the amsmath package for nicer and easier math
typesetting; have a look at its documentation (the Kopka / Daly book
also has a section about it).

* hyperref.cfg
This file contains the setup for the hyperref package which produces
a nicely hyperlinked PDF file (you could also include the hypersetup
in the root file, but that makes it kind of messy). Have a look at
the documentation of the hyperref package. You definitely have to
adjust the driver name (dvipdfm, pdftex, dvips, etc.) to the method
you want to use to produce the PDF. I found it useful to use
"hypertex" while editing the document with MiKTeX and DVI preview
(in which case there are working links in the DVI file when viewed
with the Yap program) and switching to "dvipdfm" for the final
version from which the PDF was produced (in which case the links
work fine in the PDF).


----------------------------------------------------
Modifications
----------------------------------------------------

Holger Fleckenstein, June 2007

I have changed the style file (usbthesis.cls) and thesis template
(thesis.tex) to fit the new guidelines for electronic submission.

You can now choose between the options "sig" and "nosig" (default) for
the documentclass: 

sig: The signature page contains the traditional lines for the orignal
page. The abstract page contains the advisors name (for separate abstract)

nosig: Complete style for electronic submission. No advisor name on
abstract, no signature lines, instead the dean's name is included as
specified in thesis.tex.

Also new: You also have to specify the type of thesis ("phd" or
"master") in the documentclass.
