# .cshrc, From the A/UX 3.0 Distribution, modified for MAPMAKER.
# The file is run everytime you invoke C-shell (e.g. from COmmandShell)
# Feel free to make modifications and add your own stuff...

# It is most important to set MAPM_LIB to name the directory where you put the
# MAPMAKER help files and executables. Likely places for this are either 
# "/users/Name/Mapmaker" or "/usr/local/bin".  Make sure you get the 
# capitalization exactly correct and that you have no spaces inside the
# quotation marks.

setenv MAPM_LIB "/users/Guest/Mapmaker"

# If MAPM_LIB is correct, then the lines below stick it in your path.  You don't need
# to do anything, you should just be able to type "mapmaker" or "qtl" and it'll work.

# The follow is (and should be kept as) one long line.
setenv PATH ".:/bin:/usr/bin:/usr/ucb:/mac/bin:/etc:/usr/etc:/usr/local/bin:/usr/bin/X11:/usr/local/Gnu:$MAPM_LIB"
rehash

# You can ignore everything below this point, if you wish.
# Print a blank line...

echo " "

# Set the shell prompt...
if ($?prompt) then				# interactive shell
	set prompt="\! A/UX $cwd% "
endif

# handy aliases
alias cd    'chdir \!*; set prompt="\! A/UX $cwd% "'
alias pushd 'pushd \!*; set prompt="\! A/UX $cwd% "'                            
alias popd  'popd  \!*; set prompt="\! A/UX $cwd% "'                            
alias dir   'ls'
alias del   'rm'
alias copy  'cp'
alias ren   'mv'
alias ls    'ls -FC'
alias ll    'ls -l'
alias h     'history'
alias print 'lpr -h \!*'
alias list  'lpr -h -p \!*'

set history = 30	# History stores last 30 commands
set notify			# provide notification of job completion

# The following aliases and settings can provide extra security
# To use, remove the leading #'s
# a	rm	rm -i		# prompt before removing any file
# a	cp	cp -i		# prompt before overwriting any file with cp
# a	mv	mv -i		# prompt before overwriting any file with mv
# set noclobber		# forbid use of > to accidentally overwrite files
# set ignoreeof		# Can't quit shell with control-D, only 'exit'
