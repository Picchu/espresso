from model_def import *
from saga_funcs import *

#################  Links ############################

class Include(Derive):
    def __repr__(self):
        return "<Include {0} -> {1}>".format(self.start, self.end)
    __mapper_args__ = {'polymorphic_identity': 'include'}


class Copy(Derive):
    def __repr__(self):
        return "<Include {0} -> {1}>".format(self.start, self.end)
    __mapper_args__ = {'polymorphic_identity': 'copy'}

############################# DATA SUBCLASSES ################################
# data class containing parameters for a Quantum Espresso calculation
# it does not include info on the cell itself, since that will be taken from
# a Struc object
class QEInput(Data):
    def __repr__(self):
        return "qeinput({})".format(self.name)
    __mapper_args__ = {'polymorphic_identity': 'qeinput'}

# data class containing information about a parsed Quantum Espresso
# output file. Currently parsing is done with espresso module for python
# compatible with ASE (https://github.com/zhongnanxu/espresso)
class QEOutput(Data):
    def __repr__(self):
        return "qeoutput({})".format(self.name)
    __mapper_args__ = {'polymorphic_identity': 'qeoutput'}

# data class containing information about a specific structure
# it allows to build a full cell for a given calculation
class Struc(Data):
    def __repr__(self):
        return "struc({})".format(self.name)
    __mapper_args__ = {'polymorphic_identity': 'struc'}

# data class which acts as a container for different structures
# representative of the same chemical composition
class Material(Data):
    def __repr__(self):
        return "material({})".format(self.name)
    __mapper_args__ = {'polymorphic_identity': 'material'}

# data class to store a MD trajectory
class Traj(Data):
    def __repr__(self):
        return "traj({})".format(self.name)
    __mapper_args__ = {'polymorphic_identity': 'traj'}

# data class used to store the path to the locations of raw files
# i.e. input/output files for Quantum Espresso
class FileData(Data):
    def __repr__(self):
        return "file({})".format(self.name)
    __mapper_args__ = {'polymorphic_identity': 'file'}

class ValueData(Data):
    def __init__(self, value):
        #self.value = value      # need to be careful of non-mapped attributes
        content = {'value': value}
        content_hash = hashlib.sha224(str(content)).hexdigest()
        super(ValueData, self).__init__(content=content, hash=content_hash )

    def __repr__(self):
        return "Value ({})".format(self.content['value'])
    __mapper_args__ = {'polymorphic_identity': 'valuedata'}

    @property
    def value(self):
        return self.content['value']

class IntData(ValueData):
    def __init__(self, value):
        if isinstance(value, int):
            super(IntData, self).__init__(value)
        else:
            raise Exception("Not integer data")
    def __repr__(self):
        return "Int ({})".format(self.content['value'])
    __mapper_args__ = {'polymorphic_identity': 'integer'}

class FloatData(ValueData):
    def __init__(self, value):
        if isinstance(value, float):
            super(FloatData, self).__init__(value)
        else:
            raise Exception("Not float data")
    def __repr__(self):
        return "Float ({})".format(self.content['value'])
    __mapper_args__ = {'polymorphic_identity': 'float'}

# TODO weak reference??
class StringData(ValueData):
    def __init__(self, value):
        if isinstance(value, str):
            super(StringData, self).__init__(value)
        else:
            raise Exception("Not string data")
    def __repr__(self):
        return "String ({})".format(self.content['value'])
    __mapper_args__ = {'polymorphic_identity': 'string'}

class DictData(ValueData, dict):
    def __init__(self, value):
        if isinstance(value, dict):
            super(DictData, self).__init__(value)
        else:
            raise Exception("Not dict data")
    def __repr__(self):
        return "Dict ({})".format(self.content['value'])
    __mapper_args__ = {'polymorphic_identity': 'dict'}

class ListData(ValueData, list):
    def __init__(self, value):
        if isinstance(value, list):
            super(ListData, self).__init__(value)
        else:
            raise Exception("Not list data")
    def __repr__(self):
        return "List ({})".format(self.content['value'])
    __mapper_args__ = {'polymorphic_identity': 'list'}

# more general Parameter class: input instances for specific
# programs should inherit from this one
class Param(Data):
    def __repr__(self):
        return "Param ({})".format(self.content)
    __mapper_args__ = {'polymorphic_identity': 'param'}

# data class which contains a set of Pseudo Potentials
class PPSubset(Data):
    def __repr__(self):
        return "ppsubset({})".format(self.name)
    __mapper_args__ = {'polymorphic_identity': 'ppsubset'}

# data class to store information about a given Pseudo Potential
class PseudoPot(Data):
    def __repr__(self):
        return "PseudoPot ({})".format(self.content)
    __mapper_args__ = {'polymorphic_identity': 'pseudopot'}

# data class to represent a given parsing scheme used to parse
# output files.
class Parser(Data):
    def __repr__(self):
        return "Parser ({})".format(self.name)
    __mapper_args__ = {'polymorphic_identity': 'parser'}

# object returned by an observer calc, i.e. a calc that
# periodically checks the output of another calc, killing
# it if some conditions are met
class Observed(Data):
    def __repr__(self):
        return "Observed ({})".format(self.name)
    __mapper_args__ = {'polymorphic_identity': 'observed'}

####################### CODE and CALC #########################

# instance of a Quantum Espresso calculation
class QECalc(Calc):
    def __repr__(self):
        return "QECalc [{}]".format(self.name)

    # test method which kills the calc
    # currently instead of killing the calc,
    # simply prints a message on screen
    # it is sufficient to replace that statement with
    # j.cancel() to actually kill the calc. Currently
    # it assumes all QE calcs are run as saga tasks,
    # both locally and remotely
    def kill_itself(self, username, host, js_type):
        # get id of the saga jobs associated to the calc
        # these are stored as part of the content field
        jobs_id = self.content.get('jobs_id')

        # retrieve job service associated to the saga task
        # currently we use one job service per task type,
        # so we don't have the problem of retrieving
        # the right job service instance
        job_service = get_js(username, host, js_type)

        # no jobs associated with this calc, nothing to do
        if jobs_id is None:
            return

        # loop over jobs and kill them
        for j_id in jobs_id:
            j = job_service.get_job(j_id)
            print "Job checked"

        return

    __mapper_args__ = {'polymorphic_identity': 'qecalc'}

# instance of a calculation executed by materials project
# useful to distinguish our own results from the ones pulled
# from materials project using pymatgen
class MPCalc(Calc):
    def __repr__(self):
        return "MPCalc [{}]".format(self.name)
    __mapper_args__ = {'polymorphic_identity': 'mpcalc'}

# generation of the convex hull used in stability analysis
class HullCalc(Calc):
    def __repr__(self):
        return "HullCalc [{}]".format(self.name)
    __mapper_args__ = {'polymorphic_identity': 'hullcalc'}

# this calc groups together a set of Data objects, i.e.
# generating a Pseudo Potential set, or a set of energy
# calculations outputs to be used as inputs for stability analysis
class CombineCalc(Calc):
    def __repr__(self):
        return "CombineCalc [{}]".format(self.name)
    __mapper_args__ = {'polymorphic_identity': 'combinecalc'}

class Wfunc(Code):
    def __repr__(self):
        return "Wfunc [{}]".format(self.name)
    __mapper_args__ = {'polymorphic_identity': 'wfunc'}

# currently associated to calcs in QE workflows (QECalc)
class QEplugin(Code):
    def __repr__(self):
        return "QE binary [{}]".format(self.name)
    __mapper_args__ = {'polymorphic_identity': 'qeplugin'}

# code associated to CombineCalc
class Combine(Code):
    def __repr__(self):
        return "Combine [{}]".format(self.name)
    __mapper_args__ = {'polymorphic_identity': 'combine'}

#class Join(Calc):
#    def ....
