== install ==

install with `pip install .` inside a venv


== Usage ==

This solver takes two CSV files see the example directory for the structure of the file.
And will try to come up with the most ideal assigment of fibers to fiber runs. There
are some assumptions made in about the runs / layouts so read the source code and change
some constants if this is also relevant for the usecase where one wants to use it for.

== Function ==

In the basis it uses or-tools with a set of constraints to make sure each fiber is assigned at
most once to a link with constraints that the fiber is long enough and has enough cores.
Ontop of that there are specific extras to make sure it does not selects a to big or bulky fiber.
