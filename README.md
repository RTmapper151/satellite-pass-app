# NAU Capstone

Creating a GIS app that pulls sat passovers in a given day for a given area.

LATEST BRANCH: data pull_exp

I built out three main scripts: `downloader.py`, which pulls and caches TLE files from CelesTrak; `analysis.py`, which handles satellite pass detection and plotting; and `app.py`, which ties everything together in a Streamlit interface.

I ran into an issue where switching satellite groups—like `resource` versus `scientific`—wasn’t changing the results. Turned out, some groups on CelesTrak either don’t exist or are mislabeled. The TLE URLs were failing quietly or returning fallback data, which made it look like everything was working when it wasn’t. I fixed this by adding a line in the app that shows exactly where the TLE file was pulled from—either freshly downloaded or loaded from cache. That check helps debug any group-related issues quickly.

I added group descriptions to the dropdown so when someone picks a satellite group, a short blurb appears to explain what those satellites are typically used for. Under the swath width and time interval sliders, I also added explanations: wider swaths cover more area but may reduce precision, and larger intervals mean fewer checks, which might miss narrow passes.

There’s also a note now explaining why the tool only looks at one day—it’s because most of the satellites we’re analyzing are in low Earth orbit (LEO), and their positions change quickly. Anything longer than a day risks being inaccurate.

Finally, I integrated a simple basemap using `contextily` so the plots now include country and city labels. It's still a static map, but it's a big improvement for visual clarity. The app saves a `.txt` file with the list of passing satellites and a `.png` image of the map, and both can be downloaded directly from the Streamlit interface. Everything is a lot more functional and user-friendly now.

