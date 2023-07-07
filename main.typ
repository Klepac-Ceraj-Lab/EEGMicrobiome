#import "template.typ": *

#show: project.with(
  title: "Early infant microbiobial metabolism alters brain development",
  authors: (
    "Kevin S. Bonham",
    "Emma Margolis",
    "Laurel Gabard-Durnam",
    "Vanja Klepac-Ceraj",
  ),
  // Insert your abstract after the colon, wrapped in brackets.
  // Example: `abstract: [This is my abstract...]`
  abstract: [
    Infancy is a critical window for brain development,
    and that development is shaped by numerous environmental factors.
    Emerging evidence suggests that the gut microbiome,
    which also undergoes extensive developmental changes in early life,
    may alter brain development through the metabolism
    of neuroactive compounds. Here, we show that ...
],
)

// We generated the example code below so you can see how
// your document will look. Go ahead and replace it with
// your own content!

= Introduction
#lorem(60)

= Results
#lorem(100)

#figure(caption: [Figure 1 - A caption])[#image(width: 50%, "placeholder.png")]

= Discussion
#lorem(50)

= Methods
#lorem(40)

= References
#lorem(100)
