# Python port of hp.c, edited for readability and with some extra features,
# most notably a graphical representation of what's happening.
# About a factor 100 slower, though, even without the graphics!
#
# By Mattias Sjö, 2022

import numpy as np
import random as rand
import operator as op
import os
import sys
import signal
import time
from timeit import default_timer as timer
from copy import copy

#>>> Values with "<--- ADJUST" can be adjusted if needed <<<#

# Number of iterations
N_ITER = 1_000_000  # <--- ADJUST 
# How many iterations between each printing of the output
PRINT_PERIOD = 1_000 # <--- ADJUST
# The chain is manipulated using two kinds of moves:
#  "pivots" where one half is pivoted relative to the other,
#  and "kinks" where a small portion of the chain has its shape randomized while the rest is kept fixed.
# Kinks will be done on portions ranging from a single bead up to the following length.
MAX_KINK_LENGTH = 2 # <--- ADJUST
# Pivots are relatively expensive to perform, but can modify the chain more dramatically.
# This can be set to make more pivot moves relative to the "default" amount, which is
#  "once after approximately one kink per bead and per kink-length".
PIVOT_AMOUNT = 1   # <--- ADJUST   


# If seed is empty, each simulation will be uniuqe.
# Enter a seed like rand.seed(42) to make it reproducible
# - the same seed will result in the same simulation.
rand.seed()  # <--- ADJUST

# If true, we will draw a snapshot of the chain about once a second 
#  (but no more often than output is written).
DRAW = True  # <--- ADJUST

# If true, we will draw the chain after every step of the simulation
#  (but at most MAX_FPS times a second)
# This is of course very slow, but gives great insight into how things work.
# Other output is disabled since you won't run a full simulation this way anyways.
STEP = True  # <--- ADJUST
MAX_FPS = 30 # <--- ADJUST

# If true, the chain is drawn using ncurses, which allows for very efficient drawing and flicker-less
#  display up and exceeding what FPS your screen can manage.
# However, it is not available on all platforms, so it should be turned off when running e.g.
#  under Windows or in a browser-based Jupyter Notebook.
# If false, a much more primitive variant of graphics-drawing is used, and your screen is likely
#  to flicker unless MAX_FPS is set quite low.
ADVANCED_GRAPHICS = True # <--- ADJUST
# Should not be enabled if we're not drawing anything
if not (DRAW or STEP):
    ADVANCED_GRAPHICS = False

# Loads of global variables - blame Anders!
# These are given proper values in init() - I just put them here for some semblance of clarity.
sequence = ''
n_beads = 0
origin = 0
temp = 0.0
beta = 0.0
bead_coords = []
occupied = []
energy = 0

# For graphics
if DRAW or STEP:
    if ADVANCED_GRAPHICS:
        import curses
        
        # Settings for curses
        stdscr = None
        H_COLOR = 1
        P_COLOR = 2
        HH_COLOR = 3

    else:
        # ANSI colour codes for graphics
        H_COLOR = '\033[31m'
        P_COLOR = '\033[34m'
        HH_COLOR = '\033[32m'
        RESET_COLOR = '\033[0m'

def main():
    """Setup and run the full simulation, including output etc."""
    
    # Ensure we exit gracefully when the simulation is interrupted,
    #  since it is quite likely that it will take too long and be interrupted by the user
    def exit_gracefully(sig, frame):
        # NOTE: the two arguments are needed for compatibility but aren't used.
        print("Simulation cancelled!\n")
        if ADVANCED_GRAPHICS:
            curses.endwin()
        sys.exit()
    signal.signal(signal.SIGINT, exit_gracefully)
        
    # Initiate all variables
    init()
    
    
    # Initialise graphics
    if ADVANCED_GRAPHICS:
        global stdscr
        stdscr = curses.initscr()
        
        curses.start_color()
        curses.init_pair(H_COLOR, curses.COLOR_RED, curses.COLOR_BLACK)
        curses.init_pair(P_COLOR, curses.COLOR_BLUE, curses.COLOR_BLACK)
        curses.init_pair(HH_COLOR, curses.COLOR_GREEN, curses.COLOR_BLACK)
        
        # Turn off the cursor
        curses.curs_set(0)
        
    if DRAW or STEP:
        # Initialise timer to limit redraw frequency
        last_drawn = timer()
        
    # Open output file and run simulation
    with open("output_{:.3}.dat".format(temp), 'w') as dat:
        print("Starting simulation... (NOTE: this may take a while!)")
        

        for it in np.arange(N_ITER):
            
            # Perform a sequence of small, cheap and easily accepted kink operations
            for length in np.arange(1, 1+MAX_KINK_LENGTH):
                for i in np.arange((n_beads+1 - length) // PIVOT_AMOUNT):
                    
                    if kink(length) and STEP :
                        # Limit simulation speed so graphics can keep up
                        elapsed = timer() - last_drawn
                        if elapsed < 1/MAX_FPS:
                            time.sleep(1/MAX_FPS - elapsed)
                        
                        last_drawn = timer()
                        draw()
                        
                
            # Do a big, expensive pivot operation
            if pivot() and STEP :
                # Limit simulation speed so graphics can keep up
                elapsed = timer() - last_drawn
                if elapsed < 1/MAX_FPS:
                    time.sleep(1/MAX_FPS - elapsed)
                        
                last_drawn = timer()
                draw()
            
            # Print output
            if not STEP and (it % PRINT_PERIOD == 0):
                if DRAW and (timer() - last_drawn) > 1:
                    draw()
                    
                print(f"{it}\t{-energy}", file=dat)
                print("{:>3} % done ({} iterations)".format( int( 100 * it/N_ITER ), it))
    
    # Reset graphics
    if ADVANCED_GRAPHICS:
        curses.endwin()

    print("Simulation terminated successfully!\n")
    
def init():
    """Take user input and initialise all global variables."""
        
    global sequence, n_beads, origin, temp, beta, bead_coords, occupied, energy
    
    print("""
----------------------------------
-*- HP model on square lattice -*-
----------------------------------
-   Authors: Anders Irbäck,      -
-            Daniel Nilsson,     -
-            Mattias Sjö,        -
-    and others lost to time.    -
----------------------------------
-   Ported to Python 2022        -
-    by Mattias Sjö              -
----------------------------------
""")
    
    sequence = input("Please enter sequence of H's and P's (defaults to the one in the lab manual): ")
    
    # Default sequence
    if not sequence:
        sequence = "HHHHHPPPPHPPHPHPHPPHPHPHPHP"
        print(f"    (defaulted to {sequence})")
    # Check for bad sequence
    elif any([b != 'H' and b != 'P' for b in sequence]):
        exit("Invalid sequence!")
        
    n_beads = len(sequence)
    
    temp = float( input("Please enter the temperature (dimensionless): ") )
    beta = 1/temp
    
    # Defines a grid that marks which positions are occupied 
    #  (carrying i+1 when bead i occupies it, so that 0 means "empty")
    #  and a list of bead coordinates.
    
    # We make the grid big enough that the chain can extend fully in all directions.
    # The first bead in the chain is placed in the centre of the grid.
    
    origin = n_beads+1
    bead_coords = [(origin+i, origin) for i in np.arange(n_beads)] 
    occupied = np.zeros((2*origin+1, 2*origin+1), dtype=np.int32)
    for i, (x,y) in enumerate(bead_coords):
        occupied[x,y] = i+1;
    
    # No HH bonds initially since the chain is straight
    energy = 0   
        
    
def neighbours(x, y):
    """Produce a list all coordinates orthogonally adjacent to (x,y)."""
    
    return [(x + dx, y + dy) for dx,dy in zip([1, 0, -1, 0], [0, 1, 0, -1])]

def HH_bonds(bead, valid=None):
    """Count the number of HH bonds formed by a bead.
    
    Arguments:
    bead    -- the index of the bead
    valid   -- the range of beads to which bonds should be counted (default: all)
                (this is used when some beads are meant to be moved soon and therefore
                shouldn't be counted).
    """
    
    if not valid:
        valid = range(n_beads)
    bonds = 0
    if sequence[bead] == 'H':
        for xy in neighbours(*bead_coords[bead]):
            neigh = occupied[xy]-1
            
            if neigh >= 0 and sequence[neigh] == 'H' and neigh != bead-1 and neigh != bead+1 and neigh in valid:
                bonds += 1

    return bonds

def random_point_around(coords, dist, avoid = None):
    """Get a random unoccupied point at a certain distance from a point.
    
    Returns None if no such point exists.
    
    The distance is "Manhattan distance", i.e. the total number of steps along the x and y axes,
    not the usual Euclidean distance. Therefore, the points will be uniformly distributed on the
    edge of a square diamond, not a circle.
    
    Arguments:
    coords  -- coordinates of the point as a (x,y) tuple.
    dist    -- the distance in steps.
    avoid   -- the coordinates of a point that should not be returned, despite being unoccupied
                (this is used when we want to randomly move a point and exclude its original location).
    """
    
    points = set()
    x,y = coords
    
    for d in range(0,dist+1):
        for dx,dy in [(d, dist-d), (-d, dist-d), (d, d-dist), (-d, d-dist)]:
            if not occupied[ x+dx, y+dy ]:
                points.add( (x+dx, y+dy) )
                
    points.discard(avoid)            
    
    if not points:
        return None
    
    return rand.sample(points, 1)[0]

def move(bead, coords) -> None:
    """Move a bead to new coordinates, updating all relevant data about it."""
    
    old_coords = bead_coords[bead]
    
    # Only mark square as unoccupied if another bead hasn't already moved there
    if occupied[old_coords] == bead+1:
        occupied[old_coords] = 0
        
    bead_coords[bead] = coords
    occupied[coords] = bead+1

def reset(start, old_coords) -> None:
    """Reset a part of the chain to a previous location.
    
    Arguments:
        start       -- the index of the first bead to be reset.
        old_coords  -- an array of coordinates to which the beads should be returned.
                        Its length is the number of consecutive beads to be reset.
    """
    
    for i, coords in enumerate(old_coords):
        move(start+i, coords)
        

def kink(length) -> bool:
    """
    Modify a small portion of the chain without disturbing the rest.
    
    This is a very cheap operation that is good for fine-tuning a chain
    that is already in a very densely folded sate, since such chains
    typically don't accept pivots.
    
    Arguments:
        length  -- the number of consecutive beads to be moved.
     
    Returns True if a move was possible and was accepted.
    """
    
    global energy
    
    # It's pointless to kink large parts of small chains.
    if length >= n_beads // 2:
        return False
    
    # Randomly select two beads delimiting a subchain with the given length
    # NOTE: we do not kink the beginning of the chain, since that would be a hassle.
    #  Question for you: does this introduce any bias?
    bead1 = rand.randint(1, n_beads - length)
    bead2 = bead1 + length - 1
    
    # End of the chain: randomly pick a new position for the endpoint, and then
    #  proceed as with a middle-of-the-chain kink.
    if bead2 == n_beads - 1:
        
        # Pick a random new position for bead 2
        x1,y1 = bead_coords[bead1-1]
        xy2 = random_point_around((x1,y1), length, avoid=bead_coords[bead2]) 
        
        # Failed to find any possible move
        if not xy2:
            return False
        
        x2,y2 = xy2
        
        # Pick the most straightforward way to go from bead 1 to bead 2
        steps = [(np.sign(x2-x1), 0)] * abs(x2-x1) + [(0, np.sign(y2-y1))] * abs(y2-y1)
        
        # Shuffle it to get a new randomized path
        rand.shuffle(steps)
    
    # Middle of the chain: keep the beads just before bead1 and after bead2 fixed, and 
    #  randomly change how the subchain forms a path between them.
    else:
        x1,y1 = bead_coords[ bead1-1 ]
        x2,y2 = bead_coords[ bead2+1 ]
        
        # No changes are possible on a straight line
        if x1 == x2 or y1 == y2:
            return False
        
        # Determine the currently taken path as an array of steps
        steps = [(b2[0] - b1[0], b2[1] - b1[1]) for b1,b2 in zip(bead_coords[bead1-1:bead2+1], bead_coords[bead1:bead2+2])]
        
        # Shuffle the steps, repeating until we get something distinct from the original
        oldsteps = copy(steps)
        while steps == oldsteps:
            rand.shuffle(steps)
            
        # Remove the last step since it just leads to the fixed bead after the moved beads.
        steps.pop()
    
    deltaE = 0
    
    # Save old coordinates of the beads being moved in case we need to reset
    old_coords = copy( bead_coords[bead1:bead2+1] )
    
    x,y = bead_coords[bead1-1]
    for (dx,dy), bead in zip(steps, range(bead1, bead1 + len(steps))):
        
        # Get new coordinates, check for collision
        x += dx
        y += dy
        if occupied[x,y]:
            # For efficiency, only reset beads that have been moved
            reset(bead1, old_coords[:bead - bead1])
            return False
        
        # Count broken bonds
        deltaE += HH_bonds(bead)
        
        move(bead, (x,y))
        
        # Count new bonds
        # NOTE: we may form bonds with beads that have yet to be moved.
        # This is not a problem, since those bonds will be broken right away.
        deltaE -= HH_bonds(bead)
            
    # Check if move was accepted according to thermodynamic probability
    if rand.random() < np.exp( -beta * deltaE ):
        energy += deltaE            
        return True
    # Otherwise reset to old coordinates
    else:
        reset(bead1, old_coords)
        return False
        
       
def pivot() -> bool:
    """
    Choose a "pivot point" and rigidly rotate/reflect part of the chain around it.
    
    This allows for drastic restructuring of the chain, but is a relatively expensive
    operation. Also, it is typically bad at moving between low-energy states since it
    tends to visit high-energy states between them.
    
    Returns True if a move was possible and was accepted.
    """
    
    global energy
    
    # Randomly choose a pivot (not the ends - that would be a trivial move)
    pivot_bead = rand.randint(1, n_beads-1)
    
    # Randomly choose a transformation
    # NOTE: there are 8 possible transformations of this kind, but one is trivial so we only pick between 7
    transf = rand.choice([
        lambda x,y: (-y, x), # rotate 90 deg
        lambda x,y: (-x,-y), # rotate 180 deg
        lambda x,y: ( y,-x), # rotate 270 deg
        lambda x,y: (-x, y), # mirror through y-axis
        lambda x,y: ( x,-y), # mirror through x-axis
        lambda x,y: ( y, x), # mirror through line y=x
        lambda x,y: (-x,-y)  # mirror through line y=-x
        ])

    # Copy the old coordinates in case we need to reset
    old_coords = copy( bead_coords[pivot_bead+1:] )
    
    px,py = bead_coords[pivot_bead]
    new_coords = [tuple(map(op.add, (px, py), transf(x-px, y-py))) for x,y in old_coords]
    
    deltaE = 0
    
    for i, coords in enumerate(new_coords):
        
        bead = pivot_bead+1 + i
        
        # Check collision, but only with beads not being moved
        # (moving beads can't collide with each other, or with the pivot)
        if occupied[coords]-1 in range(0,pivot_bead):
            # NOTE: must reset all beads, since some may have been overwritten by another bead
            #  despite not being moved itself.
            reset(pivot_bead+1, old_coords)
            return False
        
        # Count broken bonds, but only with beads not being moved
        # (no bonds are changed within the moving part)
        deltaE += HH_bonds(bead, valid = range(pivot_bead+1))
        
        move(bead, coords)
        
        # Count new bonds, but only with beads not being moved
        deltaE -= HH_bonds(bead, valid = range(pivot_bead+1))
        
    # Check if move is accepted
    if rand.random() < np.exp( -beta * deltaE ):
        energy += deltaE 
        return True
    # Otherwise reset to old coordinates
    else:
        reset(pivot_bead+1, old_coords)
        return False
    
#-- Here begins graphics-drawing stuff, which is separate from the physics simulation. --#
       
    
def draw() -> None:
    """Draw a visual representation of the chain to the screen.
    
    Each bead is represented as the letter H or P in red and blue, respectively.
    The links in the chain are represented by white ---'s and |'s, and the HH bonds by green *'s.
    A display showing the current energy is drawn at the top of the screen.
    
    If the chain doesn't fit on the screen, it will be cropped. Cropped parts will still show up
    as a "ghost image" in grey around the edges. The ghost image on each side consists of the entire chain
    that extends off that side, projected down on top of itself. The letter 'X' is used where two or more
    different characters overlap because of that.
    
    This uses curses for efficient terminal-based graphics, and will occupy the entire
    terminal screen while running.
    """
    
    ### Define a lot of auxiliary functions here to keep them out of the global namespace
    ### Not particularly elegant, but this is not the elegant part of the code.
    
    def get_offsets(xlo,ylo, xhi,yhi) -> (int,int):
        """Compute x and y offsets that places the chain on the screen in the best possible way.
        
        This is done as follows:
            - We try to put the middle of the chain at the center of the screen if possible.
            - If this is not possible without cropping the bead, we offset it horizontally and 
            vertically as little as possible until there is no cropping.
            - If cropping would result no matter what, we revert to centering on the middle bead.
            
        The x and y offsets are computed independently, which is useful when the bead fits in
        one dimension but not the other.
        
        Arguments:
            xlo,xhi -- with the offsets set to zero, the screen would be able to fit beads with
                        x in range(xlo,xhi).
            ylo,yhi -- like xlo,xhi but for y.
        """
        
        # Find current bounds of the chain
        minmax = lambda i: (min(i), max(i))
        xmin, xmax = minmax([x for x,_ in bead_coords])
        ymin, ymax = minmax([y for _,y in bead_coords])
        
        # Find location of the middle bead and the offset that would place it at centre-screen
        x0, y0 = bead_coords[n_beads//2]
        x0offs = x0 - (xhi - xlo)//2
        y0offs = y0 - (yhi - ylo)//2
        
        # Determine the offsets as described in the docstring.
        # Note that there are two identical copies of the code here, one for x and one for y.
        
        if (xhi - xlo) < (xmax - xmin) or ((xmin - x0offs) >= xlo and (xmax - x0offs) < xhi):
            xoffs = x0offs
        elif xmin - x0offs < xlo:
            xoffs = xmin - xlo
        elif xmax - x0offs >= xhi:
            xoffs = xmax - (xhi-1)
        else:
            assert False, "There is some mistake in the cropping algorithm!"
        
        if (yhi - ylo) < (ymax - ymin) or ((ymin - y0offs) >= ylo and (ymax - y0offs) < yhi):
            yoffs = y0offs
        elif ymin - y0offs < ylo:
            yoffs = ymin - ylo
        elif ymax - y0offs >= yhi:
            yoffs = ymax - (yhi-1)
        else:
            assert False, "There is some mistake in the cropping algorithm!"
        
        return (xoffs, yoffs)
    
    ### More auxiliaries and main drawing function for advanced graphics
    if ADVANCED_GRAPHICS:
    
        def draw_energy():
            """ Draw the energy bar on the top of the screen, and return a window representing the remaining screen.
            
            The energy bar consists of a frame with the label "energy = E" where E is the current energy,
            followed by a bar graph-style display of one green * for each HH bond present in the chain at the moment.
            """
            
            # Find out the width available to us
            # We assume the window height is sufficient, otherwise there ain't much we can do
            scr_h, scr_w = stdscr.getmaxyx() 
            
            # Width of energy bar is screen width minus margins
            width = scr_w - 4
            # 2 lines for margins, plus whatever is needed to draw the bar
            height = 2 + 1+(-energy)//width
            
            # Window for the energy and window for the rest
            energy_win = curses.newwin(height, scr_w)
            chain_win = curses.newwin(scr_h - height, scr_w, height, 0)
            
            # Draw good-looking border and write numeric energy value on it
            energy_win.border()
            energy_win.addstr(0, 2, f" energy = {energy} ", curses.A_REVERSE)
            
            # Draw a green * for each one that appears in the chain, thus visualising the number of HH bonds
            for line in range(0, (-energy)//width):
                energy_win.addstr(1+line, 3, '*' * width, curses.color_pair(HH_COLOR))
            energy_win.addstr(1+(-energy)//width, 3, '*' * ((-energy)%width), curses.color_pair(HH_COLOR))
            
            # We are done - update the energy display and return the remaining screen
            energy_win.refresh()
            return chain_win


        def screen_coords(x, y, xoffs, yoffs) -> (int,int):
            """ Convert grid x,y coordinates to coordinates on the screen. """
            
            return ((y - yoffs)*2, (x - xoffs)*4)
        def grid_coords(l, c, xoffs, yoffs) -> (int,int):
            """ Convert screen (line,column) coordinates to grid coordinates. """
            
            return (c//4 + xoffs, l//2 + yoffs)

        def draw_bond(xy1, xy2, horiz, vert, attr=0):
            """Draw the bond between two beads. 
            
            Arguments:
                xy1     -- x,y coordinates of one bead
                xy2     -- x,y coordinates of the other bead
                horiz   -- the horizontal version of the bond, e.g. "---"
                vert    -- the vertical version of the bond, e.g. "|"
                attr    -- the curses attributes with which to draw the bond
            """
            
            nonlocal lmax,cmax, xoffs,yoffs, window
            
            l1,c1 = screen_coords(*xy1, xoffs,yoffs)
            l2,c2 = screen_coords(*xy2, xoffs,yoffs)
            
            # We draw things character-by-character to make things easier for draw_element
            if l1 == l2:
                for i,ch in enumerate(horiz):
                    draw_element(l1, min(c1,c2) + 1 + i, ch, attr)
            elif c1 == c2:
                for i,ch in enumerate(vert):
                    draw_element(min(l1,l2) + 1, c1 + i, ch, attr)
            else:
                assert False, "There is some mistake in the drawing algorithm!"
                
        def draw_element(l, c, ch, attr=0):
            """Draw a character on the screen.
            
            If the given coordinates are off-screen, it will be delegated to draw_ghost appropriately.
            
            Arguments:
                l,c  -- the screen coordinates
                ch   -- the character
                attr -- the curses attributes
            """
            
            nonlocal lmax,cmax, window
            
            # Test if we're on-screen
            lgood = (l in range(1, lmax-1))
            cgood = (c in range(1, cmax-1))
            
            if lgood and cgood:
                window.addstr(l,c, ch, attr)
                
            elif lgood or (not cgood and (l == 0 or l == lmax-1)):  
                draw_ghost(l, (0 if c < 1 else cmax-1), ch)
                    
            elif cgood or (c == 0 or c == cmax-1):
                draw_ghost((0 if l < 1 else lmax-1), c, ch)
                    
                    
        
        def draw_ghost(l,c, ch):
            """Draw a ghost image as explained above."""
            
            if ch == ' ':
                pass
            
            nonlocal window
            
            # NOTE: this is not really a char, but an int. Its lowest 8 bits correspond to a char.
            oldch = window.inch(l,c) & 0b11111111

            # curses gives an error after writing to the lower right corner for some reason.
            # This try-except is here to ignore that.
            try:
                if oldch == ord(' '):
                    window.addstr(l,c, ch, curses.A_DIM)
                elif oldch != ord(ch):
                    window.addstr(l,c, 'X', curses.A_DIM)
            except curses.error:
                pass
            
        ### Now, let's proceed with the main function.
            
        # Erase previous output. 
        # curses should do this quite efficiently, so there is no need for us to figure out what actually needs
        # to be erased and what could be reused.
        stdscr.erase()
        
        # Draw a small window showing the energy, and return the remainder for the chain to be drawn on
        window = draw_energy()
        
        # Get size of window
        # To avoid confusion, we'll use (l,c) (i.e. line,column) for screen coordinates and (x,y) for chain coordinates
        lmax, cmax = window.getmaxyx()
        # Fit the chain on the screen (leaving a margin of 1 for "ghost images"
        xoffs, yoffs = get_offsets( *grid_coords(1, 1, 0, 0), *grid_coords(lmax-1, cmax-1, 0, 0) )
        
        
        # Draw all beads and their bonds
        for bead, coords in enumerate(bead_coords):
            l, c = screen_coords( *coords, xoffs, yoffs )
        
            if sequence[bead] == 'H':
                draw_element(l,c, "H", curses.color_pair(H_COLOR))
            else:
                draw_element(l,c, "P", curses.color_pair(P_COLOR))
                
            for n_coords in neighbours(*coords):
                neigh = occupied[n_coords]-1
                
                # Ignore empty squares and avoid double-counting
                if neigh < bead:
                    continue
                
                if neigh == bead+1:
                    draw_bond(coords, n_coords, "---", "|")
                        
                elif sequence[bead] == 'H' and sequence[neigh] == 'H':
                    draw_bond(coords, n_coords, " * ", "*", curses.color_pair(HH_COLOR))
            
        # Update the screen
        window.refresh()

    ### Simple graphics
    else:
        
        # Clear the screen. Should be portable unless you have a weird terminal.
        print('\033c')
        
        # Print a basic energy bar
        print(f"E = {energy}\nHH:{HH_COLOR}{'*'*-energy}{RESET_COLOR}")
        
        size = os.get_terminal_size()
        # Need 4 columns for each square
        xmax = (size.columns - 2)//4  
        # Need 2 rows for each square, also compensating for the size of the energy bar
        ymax = (size.lines - (2 + (-energy)//size.columns))//2 - 1 
        # Compute offsets to fit chain, just like with advanced graphics
        xoffs, yoffs = get_offsets(0, 0, xmax, ymax)
        
        xmin = xoffs
        xmax += xoffs
        
        ymin = yoffs
        ymax += yoffs
        
        #Draw the visible part of the grid
        for y in range(ymin, ymax): 
            
            # First draw a line of links and HH bonds (except for the first row)
            if y != ymin:
                for x in range(xmin, xmax):
                    
                    # Find vertically adjacent beads
                    b1, b2 = occupied[x,y]-1, occupied[x,y-1]-1
                    if b1 >= 0 and b2 >= 0:
                        # Neighbours in chain
                        if b1 == b2+1 or b2 == b1+1:
                            print('  | ', end='')
                        # HH bond
                        elif sequence[b1] == 'H' and sequence[b2] == 'H':
                            print(f'  {HH_COLOR}*{RESET_COLOR} ', end='')
                        # No relation
                        else:
                            print('    ', end='')
                    # No neighbour at all
                    else:
                        print('    ', end='')
            
            # New line
            print('\n  ', end='', flush=True)
            
            # Now draw a line containing beads
            for x in range(xmin, xmax):
                
                b1 = occupied[x,y]-1
                if b1 >= 0:
                    # Draw bead itself
                    if sequence[b1] == 'H':
                        print(f'{H_COLOR}H{RESET_COLOR}', end='')
                    else:
                        print(f'{P_COLOR}P{RESET_COLOR}', end='')
                        
                    
                    # Draw links and bonds
                    if x != xmax-1:
                        # Find horizontally adjacent beads
                        b2 = occupied[x+1,y]-1
                        if b2 >= 0:
                            # Neighbours in chain
                            if b1 == b2+1 or b2 == b1+1:
                                print('---', end='')
                            # HH bond
                            elif sequence[b1] == 'H' and sequence[b2] == 'H':
                                print(f' {HH_COLOR}*{RESET_COLOR} ', end='')
                            # No relation
                            else:
                                print('   ', end='')
                        # No neighbour
                        else:
                            print('   ', end='')
                # No bead
                else:
                    print('    ', end='')
                    
            # New line
            print('\n', end='', flush=True)
            
        # Final endline
        print('', flush=True)

if __name__ == "__main__":
    
    if ADVANCED_GRAPHICS:
        
        # This ensures curses doesn't mess up the terminal in case of a crash.
        # We do not use the curses wrapper since it does more than we want.
        try:                
            main()
        finally:
             # Un-initialise graphics
            curses.endwin()
    
    else:
        # Normal running without curses graphics
        main()
            
