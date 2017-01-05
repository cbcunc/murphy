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
#
# attributehandler.py is a helper file for parsing the attribute field of GFF and GTF files
#
# Routines return a dictionary containing the attribute fields and
# their values
#
# version history
# 2014-10-20  S. Stiegelmeyer  convert keys to lowercase
# 2015-07-23  S. Stiegelmeyer  add support for gtf files
# 2016-10-19  S. Stiegelmeyer  rename parseAttributes to parseAttributesGFF

def parseAttributesGFF(attrib):
    alist = attrib.split(';')
    septuple=[]
    for i in range(len(alist)):
        item = alist[i].strip().split('=')
        item[0]=item[0].lower()
        if len(item) == 2:
            septuple.append(tuple(item))
    adict = dict(septuple)
    return adict

def parseAttributesGTF(attrib):
    alist = attrib.split(';')
    septuple=[]
    for i in range(len(alist)):
        item = alist[i].strip().split(' "')
        if len(item) == 2:
            item[0]=item[0].lower()
            item[1]=item[1].strip('"')
            septuple.append(tuple(item))
    adict = dict(septuple)
    return adict

def attributesfromdictGTF(adict):
    # create an attribute field from adict
    values = []
    for item in adict:
        values.append("%s \"%s\";" % (item, adict[item]))
    return ' '.join(values)
