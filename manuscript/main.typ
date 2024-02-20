#import "template.typ": *

// Take a look at the file `template.typ` in the file panel
// to customize this template and discover how it works.
#show: project.with(
  title: "Early infant microbiobial metabolism alters brain development",
  authors: (
    "Kevin S. Bonham",
    "Emma Margolis",
    "Laurel Gabard-Durnam",
    "Vanja Klepac-Ceraj",
  ),
  abstract: [
    Infancy is a critical window for brain development,
    and that development is shaped by numerous environmental factors.
    Emerging evidence suggests that the gut microbiome,
    which also undergoes extensive developmental changes in early life,
    may alter brain development through the metabolism
    of neuroactive compounds. Here, we show that ...
],
)

// Leave comments like this

= Introduction

The brain and the gut are ...

= Results

Here are some results.

#figure(caption: [Figure 1 - A caption])[#image(width: 50%, "placeholder.png")]

And some more...

= Discussion
#lorem(50)

= Methods
#lorem(40)

= References
#lorem(100)
