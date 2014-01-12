'''build an extension module
'''

module_head = '''#define PY_SSIZE_T_CLEAN
#include "Python.h"
%s

'''

module_exception ='''static PyObject *%sError;

'''

method_obj = '''static PyObject *
%s_%s(%s)
%s
'''

method_doc = '''PyDoc_STRVAR(%s_%s_doc,
%s);

'''

method_def = '''
    {"%s", %s_%s, %s, %s_%s_doc},'''

PyMethodDef = '''
static PyMethodDef %s_methods[] = {%s
    {NULL,		NULL}		/* sentinel */
};

'''

struct_seq_field = '''
'''

struct_seq_desc = '''
'''

struct_seq_type = '''
static PyTypeObject %sType;
'''

PyStructSeqField = '''
static PyStructSequence_Field %s_type_fields[] = {%s
};
'''

PyStructSeqDesc = '''
static PyStructSequence_Desc %s_type_desc = {
    "%s.%s",
    %s,
    %s_type_fields,
    %d,
};
'''

struct_seq_init_type = '''            PyStructSequence_InitType(&%sType, &%s_type_desc);
'''

struct_seq_init_add = '''        Py_INCREF(&%sType);
        PyModule_AddObject(m, "%s", (PyObject*) &%sType);
'''

PyStructSeqTail = '''
        if (!initialized) {
%s
        }
%s
        initialized = 1;
'''

module_tail = '''
#if PY_VERSION_HEX >= 0x03000000
static struct PyModuleDef %smodule = {
	PyModuleDef_HEAD_INIT,
	"%s",
	module_doc,
	-1,
	%s_methods,
	NULL,
	NULL,
	NULL,
	NULL
};

PyMODINIT_FUNC
PyInit_%s(void)
{
	PyObject *m;
	m = PyModule_Create(&%smodule);
	/*import_array();*/
	if (m == NULL)
            return NULL;
#else
PyMODINIT_FUNC
init%s(void)
{
	PyObject *m;
	m = Py_InitModule3("%s", %s_methods, module_doc);
	/*import_array();*/
	if (m == NULL)
            goto finally;
#endif
        %sError = PyErr_NewException("%s.error", NULL, NULL);
        Py_INCREF(%sError);
        PyModule_AddObject(m, "error", %sError);
'''

class Builder(object):
    def __init__(self, module_name, include_files,
                 module_localprototype, module_localfunction,
                 module_add, module_methods,
                 module_doc, module_filename,
                 module_struct_seq, use_numpy=False):
        ## if use_numpy is True, should add "numpy/npy_3kcompat.h" and "numpy/arrayobject.h" in include files
        self.use_numpy = use_numpy
        self.module_name = module_name
        self.module_doc = module_doc
        self.module_filename = module_filename
        self.include_files = ['#include %s'%i for i in include_files]
        if self.use_numpy:
            self.include_files.extend(['#include "numpy/npy_3kcompat.h"',
                                       '#include "numpy/arrayobject.h"'])
        self.include_files = '\n'.join(self.include_files)
        self.localprototype = module_localprototype
        self.localfunction = module_localfunction
        self.module_add = module_add
        self.methods = module_methods
        self.module_method = ''
        self.method_list = ''
        self.m_t = {'METH_NOARGS':'PyObject *self',
                    'METH_O':'PyObject *self, PyObject *arg',
                    'METH_VARARGS':'PyObject *self, PyObject *args',
                    'METH_VARARGS|METH_KEYWORDS':'PyObject *self, PyObject *args, PyObject *kwds'
                    }
        self.tail_struct_seq = ''
        if module_struct_seq:
            self.localprototype = 'static int initialized;\n' + self.localprototype
            init_type = ''
            init_add = ''
            for s in module_struct_seq:
                s_name = s[0]
                s_type = '%s'%(s_name.capitalize())
                self.localprototype += struct_seq_type%s_type
                f = ''
                for i in s[1]:
                    f += '''
    {"%s", "%s"},'''%i
                f += '''
    {0}'''
                self.localfunction += PyStructSeqField%(s_name, f)
                self.localfunction += PyStructSeqDesc%(s_name, module_name, s_name, s[2], s_name, len(s[1]))
                init_type += struct_seq_init_type%(s_type, s_name)
                init_add += struct_seq_init_add%(s_type, s_name, s_type)
            self.tail_struct_seq += PyStructSeqTail%(init_type, init_add)


    def head(self):      
        return module_head%self.include_files
    
    def exception(self):
        return module_exception%self.module_name
       
    def method_def(self):
        for (m_name, m_type, m_doc, m_block) in self.methods:
            self.module_method += method_obj%(self.module_name,
                                              m_name,
                                              self.m_t[m_type],
                                              m_block)
            self.module_method += method_doc%(self.module_name,
                                              m_name,
                                              m_doc)
            if (m_type == 'METH_VARARGS|METH_KEYWORDS') or (m_type == 'METH_NOARGS'):
                self.method_list += method_def%(m_name,
                                                '(PyCFunction)'+self.module_name,
                                                m_name,
                                                m_type,
                                                self.module_name,
                                                m_name)
            else:
                self.method_list += method_def%(m_name,
                                                self.module_name,
                                                m_name,
                                                m_type,
                                                self.module_name,
                                                m_name)
            
        return self.module_method+PyMethodDef%(self.module_name,
                                               self.method_list)
    
    def tail(self):
        m = '''PyDoc_STRVAR(module_doc,
"%s");
'''%self.module_doc
        n12 = (self.module_name,)*12
        m += module_tail%n12
        for a in self.module_add:
            m += '''        PyModule_AddObject(m, "%s", %s);\n'''%a
        m += self.tail_struct_seq
        m += '''
#if PY_VERSION_HEX >= 0x03000000
        return m;
#else
        finally:
        return;
#endif
}
'''
        if self.use_numpy:
            m = m.replace('/*import_array();*/', 'import_array();')
        return m

    def create(self):
        with open(self.module_filename, 'w') as w:
            w.write(self.head())
            w.write(self.exception())
            w.write(self.localprototype)
            w.write(self.localfunction)
            w.write(self.method_def())
            w.write(self.tail())
