// The project function defines how your document looks.
// It takes your content and some metadata and formats it.
// Go ahead and customize it to your liking!
#let project(title: "", abstract: [], authors: (), body) = {
  // Set the document's basic properties.
  set document(author: authors, title: title)

  // Title row.
  align(center)[
    #block(text(weight: 700, 1.75em, title))
  ]

  // Author information.
  pad(
    top: 0.5em,
    bottom: 0.5em,
    x: 2em,
    grid(
      columns: (1fr,) * calc.min(3, authors.len()),
      gutter: 1em,
      ..authors.map(author => align(center, strong(author))),
    ),
  )

  // Main body.
  set par(justify: true)

  heading(outlined: false, numbering: none, text(0.85em)[Abstract])
  abstract

  body
}
