# Copyright 2016 Suzy M. Stiegelmeyer
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0

#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

# arghandler.py is a python utility script that will support command line
# parameter handling. Given a list of sys.argv arguments a dictionary will
# be returned where the key is the agument name and the value will be
# argument value.  '-' are removed from the key name in the dictionary.
#
# For arguments that do not have a value, None is returned as the argument
# value in the dictionary.
#
# Basic error handling is implemented to check for duplicate command
# areguments and values missing a command argument, i.e,
#
#     -i tmp.txt tmp2.txt -o out.txt.
#
# tmp2.txt would be flagged as an error.
#
# If an error is detected a boolean is returned to the calling program so
# that it can handle the error.

# Version history:
# 2016-11-16    S. Stiegelmeyer Initial version
# 2017-01-05    S. Stiegelmeyer Added some comments
# 2017-07-15    C. Calloway Flake8

import sys


def argHandler(args):
    i = 1  # args[0] is the script name
    argdict = {}
    mykey = None
    defval = None
    err = False
    while i < len(args):
        # check for an argument flag
        if args[i].startswith('-') and mykey is None:
            mykey = args[i].strip('-')
        elif mykey is not None:
            if not args[i].startswith('-') and mykey not in argdict:
                argdict[mykey] = args[i].strip('-')
            elif mykey not in argdict:
                argdict[mykey] = defval
            else:
                sys.stderr.write("%s is given multiple times.\n" % (args[i]))
                err = True
                break
            mykey = None
        else:
            sys.stderr.write("%s given without command argument.\n" %
                             (args[i]))
            err = True
            break
        i += 1

    return argdict, err
