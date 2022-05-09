# Parameter choice from normal location-shift experiment

I tried parameters for choosing epsilon-maximizers of the form
a_n = a * sqrt(log(log(n))) / sqrt(n)
and parameters for contact set estimation of the form
b_n = b * log(log(n)) / sqrt(n).

The second choice matches Linton, Song & Whang (2010), and the first choice is
just a little tighter than the contact set one.

Small a was clearly better, although choosing b was more difficult.  From
looking at results when a = 0.5 only, I chose b = 3.5 since that seemed to work
better than every choice but b = 5 (sometimes) and was not at the edge of yhe
parameter space.  It also comes close to matching the constants 3-4 used by
LSW 2010.


