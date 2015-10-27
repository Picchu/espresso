#from __future__ import absolute_import

from celery_app import my_app
from model_def import *
from wftest import *
from wflib import *
from subtypes import *

# use name explicitly to ensure consistency
# with the function loaded by the worker
@my_app.task(name='lib.tasks_app.func0_task')
def func0_task(uuid, **inp):
    # in principle the caller can be None
    if uuid is not None:
        # convert from UUID to actual db object
        self = extract_entry_from_db(uuid, 'calc')
        if ((self == -2) or (self == -1)):
            print "DB entry not found!"
            return 1
    else:
        self = None

    # do the same with the DB input objects
    data_inp = {}
    for key, value in inp.iteritems():
        data_inp[key] = extract_entry_from_db(value, 'data')
        if ((data_inp[key] == -2) or (data_inp[key] == -1)):
            print "DB entry not found!"
            return 1

    # now call the work function in its first-level wrapped state
    # otherwise we would get an infinite loop
    out = func0._original(self, **data_inp)

    session = Session()

    out_copy = {}
    for key, value in out.iteritems():
        out_copy[key] = session.merge(value)

    # return the UUIDs of the outputs in dict format
    # conversion back to DB objects will be taken care of
    # in the w4 method
    # we only need the UUIDs for the keyword arguments
    uuids = extract_uuid_from_db(**out_copy)[1]

    session.close()

    return uuids
