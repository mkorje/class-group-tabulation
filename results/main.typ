#import "@preview/lilaq:0.6.0" as lq

#set page(
  width: auto,
  height: auto,
  fill: none,
  margin: (x: 4.3em, top: 0.4em, bottom: 2.6em),
)

#set align(center)

#set text(font: "Latin Modern Roman 12", size: 12pt)
#show math.equation: set text(font: "Latin Modern Math", weight: "regular")

#show: lq.cond-set(lq.grid.with(kind: "x"), stroke: none)

#show lq.selector(lq.legend): set grid(
  row-gutter: 0em,
  column-gutter: 0em,
  columns: 4,
)
#show lq.selector(lq.legend): set text(size: 10pt)

#show: lq.set-diagram(
  width: 6cm,
  xlim: (0, 512),
  xlabel: none,
  title: none,
  ylim: (auto, 1),
  xaxis: (
    mirror: none,
    subticks: none,
    ticks: (0, 127, 255, 383, 511).map(k => (k, none)),
  ),
  yaxis: (
    mirror: none,
    exponent: none,
    position: left,
  ),
  margin: 0%,
  bounds: "data-area",
  legend: none,
)

#set grid(
  row-gutter: 1em,
  column-gutter: 1.5em,
)

#show grid.cell.where(x: 1): it => {
  show: lq.set-diagram(
    yaxis: (
      position: right,
    ),
  )
  it
}

#show grid.cell.where(y: 1): it => {
  show: lq.set-diagram(
    xaxis: (
      ticks: (0, 127, 255, 383, 511)
        .map(k => (k, $#(k + 1)$))
        .enumerate()
        .map(i => {
          if i.first() == 4 {
            let (k, eq) = i.last()
            (
              k,
              [
                #place(dy: 1.2em, dx: 0.2em, $dot 2^28$)
                $512$
              ],
            )
          } else {
            i.last()
          }
        }),
    ),
    xlabel: $x$,
  )
  it
}

#show grid.cell.where(y: 0, x: 0): it => {
  show: lq.set-diagram(
    legend: (position: bottom + right, pad: 4pt),
  )
  it
}

#let primes = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31)
#let data = (
  primes
    .map(p => (str(p), csv("data/ell" + str(p) + ".csv", row-type: dictionary)))
    .to-dict()
)

#let eta(l, k: 100) = if k > 0 {
  range(1, k + 1).map(i => 1 - calc.pow(l, -i)).product()
} else {
  1
}

#let diagram(
  label,
  kind,
) = {
  let xs = range(0, 512)
  let plot(l, total, N) = lq.plot(
    xs,
    data
      .at(str(l))
      .map(x => (
        (total(x) / int(x.at(kind)))
          / (eta(l) * (calc.pow(l, -N + 1) / (l - 1)))
      )),
    mark: none,
    smooth: true,
    label: $ell = #l$,
    stroke: 1pt,
    clip: false,
  )

  page[
    #show: lq.layout

    #grid(
      columns: 2,
      ..range(1, 5).map(N => {
        let total(x) = range(N, 10)
          .map(i => int(x.at(
            kind + "_cyclic_" + str(i) + "_rank_same",
            default: 0,
          )))
          .sum()

        let (cycle, plots) = if kind == "ramified" {
          (
            (white,) + lq.color.map.petroff10.slice(1),
            (lq.plot((), (), label: []),)
              + primes.slice(1).map(l => plot(l, total, N)),
          )
        } else {
          (lq.color.map.petroff10, primes.map(l => plot(l, total, N)))
        }

        let start = if kind == "ramified" { 1 } else { 0 }
        lq.diagram(
          ylabel: label($ell$, $#N$),
          ..plots,
          cycle: cycle,
        )
      })
    )
  ]
}

#diagram(
  (l, N) => $i_(#l, #N)(x)$,
  "inert",
)

#diagram(
  (l, N) => $r_(#l, #N)(x)$,
  "ramified",
)
