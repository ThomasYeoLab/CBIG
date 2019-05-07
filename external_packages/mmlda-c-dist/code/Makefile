# (C) Copyright 2018, Nanbo Sun (nanbosun@u.nus.edu)

# This file is part of MMLDA-C.

# MMLDA-C is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your
# option) any later version.

# MMLDA-C is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA

.SUFFIXES: .c .u
CC= gcc
CFLAGS= -O3 -Wall -g
LDFLAGS= -lm

LOBJECTS= mmlda-data.o mmlda-estimate.o mmlda-model.o mmlda-inference.o utils.o cokus.o mmlda-alpha.o

LSOURCE= mmlda-data.c mmlda-estimate.c mmlda-model.c mmlda-inference.c utils.c cokus.c mmlda-alpha.c

MMLDA:	$(LOBJECTS)
	$(CC) $(CFLAGS) $(LOBJECTS) -o MMLDA $(LDFLAGS)

clean:
	-rm -f *.o
