// CPython adapters for matching, exact pool arithmetic and record validation.
// Prepared policy executes the same reference branches on numeric evidence.
// Buffer ownership and per-query vectors keep independent calls safe.
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <new>
#include <memory>
#include <vector>

extern "C" int lg_search_plan(
    const double*, const double*, std::uint32_t, const double*, const std::uint8_t*,
    const std::uint32_t*, const std::uint8_t*, const std::uint32_t*, std::uint32_t,
    const std::uint32_t*, std::uint32_t, const std::uint32_t*, const std::uint64_t*,
    const std::uint8_t*, std::uint32_t, std::uint32_t, double, std::uint8_t,
    std::uint8_t, std::uint8_t*, std::int32_t*) noexcept;

namespace {
struct Owned {
    PyObject* value;
    explicit Owned(PyObject* p) : value(p) {}
    ~Owned() { Py_XDECREF(value); }
};
constexpr const char* QUERY_CAPSULE = "lipidgate.query_frame.v1";
struct QueryFrame {
    PyObject *peaks{}, *base{};
    std::vector<double> masses, relative, raw;
    std::vector<PyObject*> raw_objects, relative_objects;
    std::vector<std::pair<std::uintptr_t,Py_ssize_t>> positions;
    double precursor{}, exclusion{}, denominator{};
    bool float_masses=true;
    ~QueryFrame() {
        for (auto* value:raw_objects) Py_DECREF(value);
        for (auto* value:relative_objects) Py_DECREF(value);
        Py_XDECREF(peaks);Py_XDECREF(base);
    }
    Py_ssize_t position(PyObject* peak) const {
        if (!positions.empty()) {
            const auto identity=reinterpret_cast<std::uintptr_t>(peak);
            auto found=std::lower_bound(positions.begin(),positions.end(),identity,
                [](const auto& item,auto key){return item.first<key;});
            return found!=positions.end() && found->first==identity ? found->second : -1;
        }
        auto** items=PySequence_Fast_ITEMS(peaks);
        auto* found=std::find(items,items+masses.size(),peak);
        return found!=items+masses.size() ? found-items : -1;
    }
};
void release_query(PyObject* capsule) {
    auto* frame=static_cast<QueryFrame*>(PyCapsule_GetPointer(capsule,QUERY_CAPSULE));
    if (frame) delete frame;else PyErr_Clear();
}
PyObject* create_query_frame(PyObject*,PyObject*,double);
constexpr const char* CAPSULE_NAME = "lipidgate.numeric_match_plan.v1";
struct Plan {
    Py_buffer columns[9]{};
    std::uint32_t records{}, programs{};
    std::vector<std::uint32_t> plan_offsets, arm_offsets;
    std::vector<std::uint64_t> masks;
    std::vector<std::uint8_t> minimums;
    ~Plan() { for (auto& b : columns) if (b.obj) PyBuffer_Release(&b); }
};
struct ModuleState {
    PyObject *peak_type{}, *float32_type{}, *float64_type{}, *match_type{}, *fragment_type{};
    PyObject *record_type{}, *spectrum_type{};
    PyObject *rule_type{}, *profile_type{};
    PyObject* sphingo_rule_type{};
    PyObject* sphingo_rule_fields[14]{};
    PyObject* match_initializer{};
    PyObject* rule_fields[16]{};
    PyObject* profile_fields[4]{};
    PyObject* record_fields[10]{};
    PyObject* fragment_fields[6]{};
    PyObject* record_descriptors[10]{};
    PyObject* fragment_descriptors[6]{};
    Py_ssize_t record_offsets[10]{}, fragment_offsets[6]{};
    PyObject* fast_names[6]{};
    PyObject* fast_descriptors[6]{};
    Py_ssize_t fast_offsets[6]{};
    mutable unsigned int fast_versions[6]{}, record_version{}, fragment_version{};
    mutable Py_ssize_t live_fast_offsets[6]{}, live_record_offsets[10]{}, live_fragment_offsets[6]{};
};
constexpr const char* RECORD_FIELDS[] = {"record_id","compound_class","lipid_name","lipid_chain_name",
    "precursor_mz","adduct","formula","polarity","fragments","metadata"};
constexpr const char* FRAGMENT_FIELDS[] = {"mz","name","fragment_type","intensity","weight","required_group"};
constexpr const char* RULE_FIELDS[] = {"lipid_class","required_groups","chain_level_min_fah",
    "allow_hg_only_if_no_fah","positive_hg_min_matches_if_no_fah","positive_hg_min_fraction",
    "positive_hg_min_matches","positive_loss_min_fraction","positive_loss_min_matches",
    "positive_hg_complete_can_resolve_chain","positive_signature_min_matches_if_no_fah",
    "allow_loss_only_if_no_fah","require_loss_with_fah_only","negative_hg_min_matches",
    "negative_require_all_fah","score_profile"};
constexpr const char* PROFILE_FIELDS[] = {"pool_weights","metric_weights","single_group_fallback_min_matches",
    "single_group_fallback_min_relative_intensity_sum"};
constexpr const char* SPHINGO_RULE_FIELDS[] = {"lipid_class","adduct","required_all","required_any_groups",
    "required_any_by_series","required_type_any_groups","required_type_count_groups",
    "required_type_count_any_groups","required_type_fraction_groups","required_name_pattern_fraction_groups",
    "only_non_hydroxy_fa","prefer_d_series_if_ambiguous","allow_hg_only_fallback","allow_fah_only_fallback"};
constexpr const char* FAST_FIELDS[] = {"mz","relative_intensity","fragment","experimental_peak","intensity","mz_error"};
unsigned int type_version(PyTypeObject* type) {
    // CPython invalidates this tag when a class or its inherited descriptors
    // change. An unassignable/zero tag always takes the complete check below.
#if PY_VERSION_HEX >= 0x030D0000
    if (!type->tp_version_tag) PyUnstable_Type_AssignVersionTag(type);
    return type->tp_version_tag;
#else
    return 0;
#endif
}
Py_ssize_t active_offset(const ModuleState* state, int index) {
    auto* type=reinterpret_cast<PyTypeObject*>(index<2 || index==4 ? state->peak_type : state->match_type);
    const auto version=type_version(type);
    if (version && state->fast_versions[index]==version) return state->live_fast_offsets[index];
    const auto offset=(type->tp_getattro==PyObject_GenericGetAttr && state->fast_offsets[index]>=0 &&
            PyDict_GetItemWithError(type->tp_dict,state->fast_names[index])==state->fast_descriptors[index])
        ? state->fast_offsets[index] : -1;
    state->fast_versions[index]=version;state->live_fast_offsets[index]=offset;
    return offset;
}
Py_ssize_t descriptor_offset(PyTypeObject* type, PyObject* descriptor) {
    if (descriptor && Py_TYPE(descriptor)==&PyMemberDescr_Type && PyDescr_TYPE(descriptor)==type) {
        auto* member=reinterpret_cast<PyMemberDescrObject*>(descriptor)->d_member;
        if (member->type==Py_T_OBJECT_EX && member->offset>=sizeof(PyObject) &&
                member->offset+sizeof(PyObject*)<=type->tp_basicsize) return member->offset;
    }
    return -1;
}
template<std::size_t N>
void active_offsets(const ModuleState* state, PyObject* model, PyObject* const (&names)[N],
                    PyObject* const (&descriptors)[N], const Py_ssize_t (&cached)[N],
                    Py_ssize_t (&result)[N]) {
    auto* type=reinterpret_cast<PyTypeObject*>(model);
    auto& seen=N==10 ? state->record_version : state->fragment_version;
    auto* live=N==10 ? state->live_record_offsets : state->live_fragment_offsets;
    const auto version=type_version(type);
    if (version && version==seen) {
        std::copy(live,live+N,result);return;
    }
    for (std::size_t i=0; i<N; ++i)
        result[i]=(type->tp_getattro==PyObject_GenericGetAttr && cached[i]>=0 &&
            PyDict_GetItemWithError(type->tp_dict,names[i])==descriptors[i]) ? cached[i] : -1;
    seen=version;std::copy(result,result+N,live);
}
PyObject* model_field(PyObject* object, PyObject* name, Py_ssize_t offset) {
    if (offset>=0) {
        auto* value=*reinterpret_cast<PyObject**>(reinterpret_cast<char*>(object)+offset);
        if (value) return Py_NewRef(value);
    }
    return PyObject_GetAttr(object,name);
}
PyObject* peak_or_match_field(PyObject* object,const ModuleState* state,int index,Py_ssize_t offset) {
    return model_field(object,state->fast_names[index],offset);
}
void release_plan(PyObject* capsule) {
    auto* plan = static_cast<Plan*>(PyCapsule_GetPointer(capsule, CAPSULE_NAME));
    if (plan) delete plan;
    else PyErr_Clear();
}
bool valid_buffer(const Py_buffer& b, char kind, Py_ssize_t width) {
    return b.ndim == 1 && b.itemsize == width && b.len % width == 0 &&
        b.format && b.format[0] == kind && b.format[1] == '\0' && PyBuffer_IsContiguous(&b, 'C');
}
PyObject* create_plan(PyObject*, PyObject* args) {
    if (PyTuple_GET_SIZE(args) != 9) {
        PyErr_SetString(PyExc_TypeError, "Expected nine validated numeric columns"); return nullptr;
    }
    auto* plan = new (std::nothrow) Plan;
    if (!plan) return PyErr_NoMemory();
    const char kinds[] = {'d','B','I','B','I','I','I','Q','B'};
    const Py_ssize_t widths[] = {8,1,4,1,4,4,4,8,1};
    for (int i = 0; i < 9; ++i) {
        if (PyObject_GetBuffer(PyTuple_GET_ITEM(args,i), &plan->columns[i], PyBUF_FORMAT | PyBUF_ND) < 0) {
            delete plan; return nullptr;
        }
        if (!valid_buffer(plan->columns[i], kinds[i], widths[i])) {
            delete plan; PyErr_SetString(PyExc_ValueError, "Invalid native match column"); return nullptr;
        }
    }
    const auto length = [&](int i) { return plan->columns[i].len / widths[i]; };
    const auto record_count = length(3), program_count = length(5) - 1;
    if (record_count > std::numeric_limits<std::uint32_t>::max() || program_count < 0 ||
        program_count > std::numeric_limits<std::uint32_t>::max() ||
        length(0) != length(1) || length(2) != record_count+1 || length(4) != record_count ||
        length(7) != length(8) || length(6) < 1) {
        delete plan; PyErr_SetString(PyExc_ValueError, "Inconsistent native match columns"); return nullptr;
    }
    auto* offsets = static_cast<const std::uint32_t*>(plan->columns[2].buf);
    auto* programs = static_cast<const std::uint32_t*>(plan->columns[4].buf);
    auto* program_offsets = static_cast<const std::uint32_t*>(plan->columns[5].buf);
    auto* arm_offsets = static_cast<const std::uint32_t*>(plan->columns[6].buf);
    bool valid = offsets[0] == 0 && offsets[record_count] == length(0) &&
        program_offsets[0] == 0 && program_offsets[program_count] == length(6)-1 &&
        arm_offsets[0] == 0 && arm_offsets[length(6)-1] == length(7);
    for (Py_ssize_t i=0; valid && i<record_count; ++i)
        valid = offsets[i] <= offsets[i+1] && programs[i] < program_count;
    for (Py_ssize_t i=0; valid && i<program_count; ++i)
        valid = program_offsets[i] <= program_offsets[i+1];
    for (Py_ssize_t i=0; valid && i<length(6)-1; ++i)
        valid = arm_offsets[i] <= arm_offsets[i+1];
    if (!valid) { delete plan; PyErr_SetString(PyExc_ValueError,"Invalid native match offsets"); return nullptr; }
    plan->records = static_cast<std::uint32_t>(record_count);
    plan->programs = static_cast<std::uint32_t>(program_count);
    // Tiny gate programs get private, immutable copies. Mutating an exported
    // Python buffer cannot turn a previously checked clause into an OOB read.
    try {
        plan->plan_offsets.assign(program_offsets, program_offsets + length(5));
        plan->arm_offsets.assign(arm_offsets, arm_offsets + length(6));
        auto* masks=static_cast<const std::uint64_t*>(plan->columns[7].buf);
        auto* minimums=static_cast<const std::uint8_t*>(plan->columns[8].buf);
        if (length(7)) plan->masks.assign(masks,masks+length(7));
        if (length(8)) plan->minimums.assign(minimums,minimums+length(8));
    } catch (const std::bad_alloc&) { delete plan;return PyErr_NoMemory(); }
    auto* result = PyCapsule_New(plan, CAPSULE_NAME, release_plan);
    if (!result) delete plan;
    return result;
}
// Return 0 for custom/unrepresentable values; callers retain the original path.
int builtin_number(PyObject* value, double& result) {
    if (PyFloat_CheckExact(value)) result = PyFloat_AS_DOUBLE(value);
    else if (PyLong_CheckExact(value)) {
        result = PyLong_AsDouble(value);
        if (PyErr_Occurred()) { PyErr_Clear(); return 0; }
        auto* roundtrip = PyLong_FromDouble(result);
        if (!roundtrip) { PyErr_Clear(); return 0; }
        const int same = PyObject_RichCompareBool(value, roundtrip, Py_EQ);
        Py_DECREF(roundtrip);
        if (same < 0) return -1;
        if (!same) return 0;
    } else return 0;
    return std::isfinite(result) ? 1 : 0;
}
PyObject* match_plan_impl(PyObject* module, PyObject* args,bool query_mode) {
    PyObject *capsule, *peaks;
    Py_ssize_t left, right;
    double tolerance,exclusion=0.;
    int use_ppm, prefilter;
    if (!(query_mode ? PyArg_ParseTuple(args,"OOnndppd",&capsule,&peaks,&left,&right,&tolerance,&use_ppm,&prefilter,&exclusion) :
                       PyArg_ParseTuple(args,"OOnndpp",&capsule,&peaks,&left,&right,&tolerance,&use_ppm,&prefilter))) return nullptr;
    auto* plan = static_cast<Plan*>(PyCapsule_GetPointer(capsule,CAPSULE_NAME));
    if (!plan) return nullptr;
    if (left < 0 || left > right || right > plan->records || !std::isfinite(tolerance) || tolerance < 0) {
        PyErr_SetString(PyExc_ValueError,"Invalid native query bounds"); return nullptr;
    }
    Owned query(query_mode ? create_query_frame(module,peaks,exclusion) : Py_NewRef(Py_None));
    if (!query.value) return nullptr;
    if (query_mode && query.value==Py_None) Py_RETURN_NONE;
    auto* frame=query_mode ? static_cast<QueryFrame*>(PyCapsule_GetPointer(query.value,QUERY_CAPSULE)) : nullptr;
    if (frame) peaks=frame->peaks;
    if (!PyList_CheckExact(peaks) && !PyTuple_CheckExact(peaks)) Py_RETURN_NONE;
    const auto count = PySequence_Fast_GET_SIZE(peaks);
    if (count > std::numeric_limits<std::int32_t>::max()) Py_RETURN_NONE;
    auto* state = static_cast<ModuleState*>(PyModule_GetState(module));
    const auto mass_offset=active_offset(state,0), intensity_offset=active_offset(state,1);
    auto** items = PySequence_Fast_ITEMS(peaks);
    try {
        std::vector<double> masses(frame ? 0 : count), intensities(frame ? 0 : count);
        bool has_float32 = false, has_builtin = false;
        for (Py_ssize_t i=0; !frame && i<count; ++i) {
            if (reinterpret_cast<PyObject*>(Py_TYPE(items[i])) != state->peak_type) Py_RETURN_NONE;
            PyObject* mass = peak_or_match_field(items[i],state,0,mass_offset);
            if (!mass) return nullptr;
            int supported = builtin_number(mass,masses[i]); Py_DECREF(mass);
            if (supported < 0) return nullptr;
            if (!supported || (i && masses[i-1] > masses[i])) Py_RETURN_NONE;
            PyObject* intensity = peak_or_match_field(items[i],state,1,intensity_offset);
            if (!intensity) return nullptr;
            PyObject* type = reinterpret_cast<PyObject*>(Py_TYPE(intensity));
            if (type == state->float32_type || type == state->float64_type) {
                has_float32 |= type == state->float32_type;
                intensities[i] = PyFloat_AsDouble(intensity);
                supported = PyErr_Occurred() ? -1 : static_cast<int>(std::isfinite(intensities[i]));
            } else {
                has_builtin = true;
                supported = builtin_number(intensity,intensities[i]);
            }
            Py_DECREF(intensity);
            if (supported < 0) return nullptr;
            if (!supported || (has_float32 && has_builtin)) Py_RETURN_NONE;
        }
        auto* offsets = static_cast<const std::uint32_t*>(plan->columns[2].buf);
        const auto base = offsets[left], stop = offsets[right];
        if (base > stop || stop > plan->columns[0].len / 8) {
            PyErr_SetString(PyExc_ValueError,"Native query offsets changed");return nullptr;
        }
        auto* programs = static_cast<const std::uint32_t*>(plan->columns[4].buf);
        for (Py_ssize_t r=left; r<right; ++r) {
            if (offsets[r] < base || offsets[r] > offsets[r+1] || offsets[r+1] > stop ||
                    programs[r] >= plan->programs) {
                PyErr_SetString(PyExc_ValueError,"Native query columns changed");return nullptr;
            }
        }
        std::vector<std::uint8_t> selected(right-left);
        std::vector<std::int32_t> matches(stop-base);
        int status;
        // Hold the GIL for this short call: exported numeric buffers remain
        // live Python arrays, so a second Python thread must not edit them
        // between the bounds checks and the kernel's reads.
        status = lg_search_plan(frame ? frame->masses.data() : masses.data(),
            frame ? frame->relative.data() : intensities.data(),static_cast<std::uint32_t>(count),
            static_cast<const double*>(plan->columns[0].buf),
            static_cast<const std::uint8_t*>(plan->columns[1].buf),offsets,
            static_cast<const std::uint8_t*>(plan->columns[3].buf),
            static_cast<const std::uint32_t*>(plan->columns[4].buf),plan->records,
            plan->plan_offsets.data(),plan->programs,plan->arm_offsets.data(),
            plan->masks.data(),plan->minimums.data(),
            static_cast<std::uint32_t>(left),static_cast<std::uint32_t>(right),tolerance,
            static_cast<std::uint8_t>(use_ppm),static_cast<std::uint8_t>(prefilter),selected.data(),matches.data());
        if (status) { PyErr_SetString(PyExc_RuntimeError,"Native query failed"); return nullptr; }
        PyObject* result = PyDict_New(); if (!result) return nullptr;
        for (Py_ssize_t r=left; r<right; ++r) {
            const auto flag = selected[r-left]; if (!flag) continue;
            PyObject* value = nullptr;
            if (flag == 2) { value = Py_None; Py_INCREF(value); }
            else {
                value = PyTuple_New(offsets[r+1]-offsets[r]);
                if (value) for (auto j=offsets[r]; j<offsets[r+1]; ++j) {
                    PyObject* hit = PyLong_FromLong(matches[j-base]);
                    if (!hit) { Py_DECREF(value); value=nullptr; break; }
                    PyTuple_SET_ITEM(value,j-offsets[r],hit);
                }
            }
            PyObject* key = value ? PyLong_FromSsize_t(r) : nullptr;
            if (!key || PyDict_SetItem(result,key,value) < 0) {
                Py_XDECREF(key);Py_XDECREF(value);Py_DECREF(result);return nullptr;
            }
            Py_DECREF(key);Py_DECREF(value);
        }
        if (query_mode) {Owned owned_result(result);return PyTuple_Pack(2,result,query.value);}
        return result;
    } catch (const std::bad_alloc&) { return PyErr_NoMemory(); }
      catch (...) { PyErr_SetString(PyExc_RuntimeError,"Native adapter failed"); return nullptr; }
}
PyObject* match_plan(PyObject* module,PyObject* args) {return match_plan_impl(module,args,false);}
PyObject* match_query_plan(PyObject* module,PyObject* args) {return match_plan_impl(module,args,true);}
int pool_number(PyObject* value, const ModuleState* state, double& result) {
    auto* type = reinterpret_cast<PyObject*>(Py_TYPE(value));
    if (!PyFloat_CheckExact(value) && !PyLong_CheckExact(value) &&
            type != state->float32_type && type != state->float64_type) return 0;
    result = PyFloat_AsDouble(value);
    return PyErr_Occurred() ? -1 : static_cast<int>(std::isfinite(result));
}

// The Python fact scope is limited to one canonical query. Check exactly the
// same scalar/peak types here without a second Python loop over the spectrum.
bool plain_spectrum(PyObject* spectrum,const ModuleState* state) {
    auto* type=Py_TYPE(spectrum);
    return reinterpret_cast<PyObject*>(type)==state->spectrum_type &&
        type->tp_getattro==PyObject_GenericGetAttr && type->tp_base==&PyBaseObject_Type &&
        !PyDict_GetItemString(type->tp_dict,"precursor_mz") && !PyDict_GetItemString(type->tp_dict,"peaks");
}

PyObject* create_query_frame(PyObject* module,PyObject* spectrum,double exclusion) {
    auto* state=static_cast<ModuleState*>(PyModule_GetState(module));
    if (!plain_spectrum(spectrum,state) || !std::isfinite(exclusion)) Py_RETURN_NONE;
    const auto mz_offset=active_offset(state,0), relative_offset=active_offset(state,1), raw_offset=active_offset(state,4);
    if (mz_offset<0 || relative_offset<0 || raw_offset<0) Py_RETURN_NONE;
    Owned precursor(PyObject_GetAttrString(spectrum,"precursor_mz"));
    Owned peaks(PyObject_GetAttrString(spectrum,"peaks"));
    if (!precursor.value || !peaks.value) return nullptr;
    if (!PyList_CheckExact(peaks.value) && !PyTuple_CheckExact(peaks.value)) Py_RETURN_NONE;
    const auto count=PySequence_Fast_GET_SIZE(peaks.value);
    if (count>std::numeric_limits<std::int32_t>::max()) Py_RETURN_NONE;
    try {
        auto frame=std::make_unique<QueryFrame>();frame->exclusion=exclusion;
        const int p=builtin_number(precursor.value,frame->precursor);
        if (p<0) return nullptr;if (!p) Py_RETURN_NONE;
        frame->peaks=PySequence_Tuple(peaks.value);if (!frame->peaks) return nullptr;
        frame->base=PyFloat_FromDouble(0.);if (!frame->base) return nullptr;
        frame->masses.resize(count);frame->relative.resize(count);frame->raw.resize(count);
        frame->raw_objects.reserve(count);frame->relative_objects.reserve(count);
        if (count>256) frame->positions.reserve(count);
        bool has_float32=false,has_builtin=false,found=false;
        for (Py_ssize_t i=0;i<count;++i) {
            auto* peak=PyTuple_GET_ITEM(frame->peaks,i);
            if (reinterpret_cast<PyObject*>(Py_TYPE(peak))!=state->peak_type) Py_RETURN_NONE;
            Owned mass(peak_or_match_field(peak,state,0,mz_offset));
            Owned relative(peak_or_match_field(peak,state,1,relative_offset));
            Owned raw(peak_or_match_field(peak,state,4,raw_offset));
            if (!mass.value || !relative.value || !raw.value) return nullptr;
            const int a=builtin_number(mass.value,frame->masses[i]);
            const int b=builtin_number(raw.value,frame->raw[i]);
            auto* relative_type=reinterpret_cast<PyObject*>(Py_TYPE(relative.value));
            const bool numpy=relative_type==state->float32_type || relative_type==state->float64_type;
            has_float32 |= relative_type==state->float32_type;has_builtin |= !numpy;
            const int c=numpy ? pool_number(relative.value,state,frame->relative[i]) : builtin_number(relative.value,frame->relative[i]);
            if (a<0 || b<0 || c<0) return nullptr;
            if (!a || !b || !c || (i && frame->masses[i-1]>frame->masses[i]) ||
                    (has_float32 && has_builtin)) Py_RETURN_NONE;
            frame->float_masses &= PyFloat_CheckExact(mass.value);
            frame->raw_objects.push_back(raw.value);Py_INCREF(raw.value);
            frame->relative_objects.push_back(relative.value);Py_INCREF(relative.value);
            if (count>256) frame->positions.emplace_back(reinterpret_cast<std::uintptr_t>(peak),i);
            if (std::abs(frame->masses[i]-frame->precursor)>exclusion && (!found || frame->raw[i]>frame->denominator)) {
                frame->denominator=frame->raw[i];Py_SETREF(frame->base,Py_NewRef(raw.value));found=true;
            }
        }
        std::sort(frame->positions.begin(),frame->positions.end());
        auto* capsule=PyCapsule_New(frame.get(),QUERY_CAPSULE,release_query);
        if (capsule) frame.release();return capsule;
    } catch (const std::bad_alloc&) {return PyErr_NoMemory();}
      catch (...) {PyErr_SetString(PyExc_RuntimeError,"Native query preparation failed");return nullptr;}
}

PyObject* query_base(PyObject*,PyObject* capsule) {
    auto* frame=static_cast<QueryFrame*>(PyCapsule_GetPointer(capsule,QUERY_CAPSULE));
    return frame ? PyTuple_Pack(1,frame->base) : nullptr;
}

PyObject* query_intensity_sums(PyObject* module,PyObject* args) {
    PyObject *matches,*capsule;
    if (!PyArg_ParseTuple(args,"OO",&matches,&capsule)) return nullptr;
    auto* frame=static_cast<QueryFrame*>(PyCapsule_GetPointer(capsule,QUERY_CAPSULE));
    if (!frame) return nullptr;
    if (!PyList_CheckExact(matches) && !PyTuple_CheckExact(matches)) Py_RETURN_NONE;
    auto* state=static_cast<ModuleState*>(PyModule_GetState(module));
    const auto peak_offset=active_offset(state,3), count=PySequence_Fast_GET_SIZE(matches);
    Owned raw(PyTuple_New(count)),relative(PyTuple_New(count));
    if (!raw.value || !relative.value) return nullptr;
    for (Py_ssize_t i=0;i<count;++i) {
        auto* match=PySequence_Fast_GET_ITEM(matches,i);
        if (reinterpret_cast<PyObject*>(Py_TYPE(match))!=state->match_type) Py_RETURN_NONE;
        Owned peak(peak_or_match_field(match,state,3,peak_offset));if (!peak.value) return nullptr;
        const auto position=frame->position(peak.value);if (position<0) Py_RETURN_NONE;
        PyTuple_SET_ITEM(raw.value,i,Py_NewRef(frame->raw_objects[position]));
        PyTuple_SET_ITEM(relative.value,i,Py_NewRef(frame->relative_objects[position]));
    }
    // Use CPython's actual sum, including 3.13 compensated float summation and
    // NumPy promotion, rather than replacing it with a C++ accumulation.
    auto* sum=PyDict_GetItemString(PyEval_GetBuiltins(),"sum");
    if (!sum || !PyCFunction_Check(sum)) Py_RETURN_NONE;
    Owned a(PyObject_CallOneArg(sum,raw.value)),b(PyObject_CallOneArg(sum,relative.value));
    if (!a.value || !b.value) return nullptr;
    return PyTuple_Pack(2,a.value,b.value);
}
PyObject* spectrum_is_builtin(PyObject* module, PyObject* spectrum) {
    auto* state=static_cast<ModuleState*>(PyModule_GetState(module));
    if (reinterpret_cast<PyObject*>(Py_TYPE(spectrum))!=state->spectrum_type) Py_RETURN_FALSE;
    const auto mass_offset=active_offset(state,0), intensity_offset=active_offset(state,4);
    if (mass_offset<0 || intensity_offset<0 || !plain_spectrum(spectrum,state)) Py_RETURN_NONE;
    Owned precursor(PyObject_GetAttrString(spectrum,"precursor_mz"));
    Owned peaks(PyObject_GetAttrString(spectrum,"peaks"));
    if (!precursor.value || !peaks.value) return nullptr;
    if ((!PyFloat_CheckExact(precursor.value) && !PyLong_CheckExact(precursor.value)) ||
            (!PyList_CheckExact(peaks.value) && !PyTuple_CheckExact(peaks.value))) Py_RETURN_FALSE;
    for (Py_ssize_t i=0;i<PySequence_Fast_GET_SIZE(peaks.value);++i) {
        auto* peak=PySequence_Fast_GET_ITEM(peaks.value,i);
        if (reinterpret_cast<PyObject*>(Py_TYPE(peak))!=state->peak_type) Py_RETURN_FALSE;
        Owned mass(peak_or_match_field(peak,state,0,mass_offset));
        Owned intensity(peak_or_match_field(peak,state,4,intensity_offset));
        if (!mass.value || !intensity.value) return nullptr;
        if ((!PyFloat_CheckExact(mass.value) && !PyLong_CheckExact(mass.value)) ||
                (!PyFloat_CheckExact(intensity.value) && !PyLong_CheckExact(intensity.value))) Py_RETURN_FALSE;
    }
    Py_RETURN_TRUE;
}

// Retain the original maximum object, including its int/float type and tie
// order. Unusual arithmetic and replaced descriptors use the Python reference.
PyObject* non_precursor_base(PyObject* module, PyObject* args) {
    PyObject* spectrum;double exclusion;
    if (!PyArg_ParseTuple(args,"Od",&spectrum,&exclusion)) return nullptr;
    auto* state=static_cast<ModuleState*>(PyModule_GetState(module));
    if (!plain_spectrum(spectrum,state)) Py_RETURN_NONE;
    Owned precursor_value(PyObject_GetAttrString(spectrum,"precursor_mz"));
    Owned peaks_value(PyObject_GetAttrString(spectrum,"peaks"));
    if (!precursor_value.value || !peaks_value.value) return nullptr;
    auto* precursor=precursor_value.value;auto* peaks=peaks_value.value;
    if ((!PyList_CheckExact(peaks) && !PyTuple_CheckExact(peaks)) || !std::isfinite(exclusion)) Py_RETURN_NONE;
    double precursor_mz;
    if (builtin_number(precursor,precursor_mz)!=1) {
        if (PyErr_Occurred()) return nullptr;Py_RETURN_NONE;
    }
    const auto mass_offset=active_offset(state,0), intensity_offset=active_offset(state,4);
    if (mass_offset<0 || intensity_offset<0) Py_RETURN_NONE;
    Owned maximum(PyFloat_FromDouble(0.));if (!maximum.value) return nullptr;
    bool found=false;
    for (Py_ssize_t i=0;i<PySequence_Fast_GET_SIZE(peaks);++i) {
        auto* peak=PySequence_Fast_GET_ITEM(peaks,i);
        if (reinterpret_cast<PyObject*>(Py_TYPE(peak))!=state->peak_type) Py_RETURN_NONE;
        Owned mass(peak_or_match_field(peak,state,0,mass_offset));if (!mass.value) return nullptr;
        double mz;const int supported=builtin_number(mass.value,mz);
        if (supported<0) return nullptr;if (!supported) Py_RETURN_NONE;
        if (std::abs(mz-precursor_mz)<=exclusion) continue;
        Owned intensity(peak_or_match_field(peak,state,4,intensity_offset));if (!intensity.value) return nullptr;
        if ((!PyFloat_CheckExact(intensity.value) && !PyLong_CheckExact(intensity.value)) ||
                (PyFloat_CheckExact(intensity.value) && !std::isfinite(PyFloat_AS_DOUBLE(intensity.value)))) Py_RETURN_NONE;
        const int greater=found ? PyObject_RichCompareBool(intensity.value,maximum.value,Py_GT) : 1;
        if (greater<0) return nullptr;
        if (greater) {Py_SETREF(maximum.value,Py_NewRef(intensity.value));found=true;}
    }
    return PyTuple_Pack(1,maximum.value);
}

PyObject* materialize_matches_impl(PyObject* module, PyObject* args,bool query_mode) {
    PyObject *peaks,*fragments,*indexes,*constructor;
    if (!PyArg_ParseTuple(args,"OOOO",&peaks,&fragments,&indexes,&constructor)) return nullptr;
    auto* frame=query_mode ? static_cast<QueryFrame*>(PyCapsule_GetPointer(peaks,QUERY_CAPSULE)) : nullptr;
    if (query_mode && !frame) return nullptr;
    if (frame) {if (!frame->float_masses) Py_RETURN_NONE;peaks=frame->peaks;}
    auto* state=static_cast<ModuleState*>(PyModule_GetState(module));
    if (constructor!=state->match_type ||
            (!PyList_CheckExact(peaks) && !PyTuple_CheckExact(peaks)) ||
            (!PyList_CheckExact(fragments) && !PyTuple_CheckExact(fragments)) ||
            (!PyList_CheckExact(indexes) && !PyTuple_CheckExact(indexes)) ||
            !state->match_initializer || !PyTuple_CheckExact(state->match_initializer) ||
            PyTuple_GET_SIZE(state->match_initializer)!=2) Py_RETURN_NONE;
    auto* type=reinterpret_cast<PyTypeObject*>(constructor);
    Owned initializer(PyObject_GetAttrString(constructor,"__init__"));if (!initializer.value) return nullptr;
    if (initializer.value!=PyTuple_GET_ITEM(state->match_initializer,0) ||
            !PyFunction_Check(initializer.value) ||
            PyFunction_GET_CODE(initializer.value)!=PyTuple_GET_ITEM(state->match_initializer,1) ||
            type->tp_new!=PyBaseObject_Type.tp_new || type->tp_alloc!=PyType_GenericAlloc || type->tp_finalize) Py_RETURN_NONE;
    const auto fragment_offset=active_offset(state,2), peak_offset=active_offset(state,3);
    const auto error_offset=active_offset(state,5), mz_offset=active_offset(state,0);
    Py_ssize_t fragment_offsets[6];
    active_offsets(state,state->fragment_type,state->fragment_fields,state->fragment_descriptors,state->fragment_offsets,fragment_offsets);
    if (fragment_offset<0 || peak_offset<0 || error_offset<0 || (!frame && mz_offset<0) || fragment_offsets[0]<0) Py_RETURN_NONE;
    const auto count=std::min(PySequence_Fast_GET_SIZE(fragments),PySequence_Fast_GET_SIZE(indexes));
    if (count>64) Py_RETURN_NONE;
    struct Inputs {
        struct Item {PyObject *fragment{},*peak{};double error{};};
        std::array<Item,64> items{};Py_ssize_t count{};
        ~Inputs() {for (Py_ssize_t i=0;i<count;++i) {Py_DECREF(items[i].fragment);Py_DECREF(items[i].peak);}}
    } inputs;
    // Validate the complete input before allocating. Unsupported public objects
    // must execute their constructor/arithmetic exactly once in the reference.
    for (Py_ssize_t i=0;i<count;++i) {
        auto* index=PySequence_Fast_GET_ITEM(indexes,i);
        if (!PyLong_CheckExact(index)) Py_RETURN_NONE;
        const auto position=PyLong_AsSsize_t(index);
        if (PyErr_Occurred()) {PyErr_Clear();Py_RETURN_NONE;}
        if (position<0) continue;
        if (position>=PySequence_Fast_GET_SIZE(peaks)) Py_RETURN_NONE;
        auto* peak=PySequence_Fast_GET_ITEM(peaks,position);
        auto* fragment=PySequence_Fast_GET_ITEM(fragments,i);
        if (reinterpret_cast<PyObject*>(Py_TYPE(peak))!=state->peak_type ||
                reinterpret_cast<PyObject*>(Py_TYPE(fragment))!=state->fragment_type) Py_RETURN_NONE;
        Owned mass(frame ? nullptr : peak_or_match_field(peak,state,0,mz_offset));
        Owned expected(model_field(fragment,state->fragment_fields[0],fragment_offsets[0]));
        if ((!frame && !mass.value) || !expected.value) return nullptr;
        if ((!frame && !PyFloat_CheckExact(mass.value)) || !PyFloat_CheckExact(expected.value)) Py_RETURN_NONE;
        auto& input=inputs.items[inputs.count++];
        input.fragment=Py_NewRef(fragment);input.peak=Py_NewRef(peak);
        input.error=std::abs((frame ? frame->masses[position] : PyFloat_AS_DOUBLE(mass.value))-PyFloat_AS_DOUBLE(expected.value));
    }
    Owned result(PyList_New(inputs.count));if (!result.value) return nullptr;
    for (Py_ssize_t i=0;i<inputs.count;++i) {
        const auto& input=inputs.items[i];
        Owned error(PyFloat_FromDouble(input.error));
        if (!error.value) return nullptr;
        auto* match=type->tp_alloc(type,0);if (!match) return nullptr;
        *reinterpret_cast<PyObject**>(reinterpret_cast<char*>(match)+fragment_offset)=Py_NewRef(input.fragment);
        *reinterpret_cast<PyObject**>(reinterpret_cast<char*>(match)+peak_offset)=Py_NewRef(input.peak);
        *reinterpret_cast<PyObject**>(reinterpret_cast<char*>(match)+error_offset)=Py_NewRef(error.value);
        PyList_SET_ITEM(result.value,i,match);
    }
    return Py_NewRef(result.value);
}
PyObject* materialize_matches(PyObject* module,PyObject* args) {return materialize_matches_impl(module,args,false);}
PyObject* materialize_query_matches(PyObject* module,PyObject* args) {return materialize_matches_impl(module,args,true);}

// Preserve the reference float conversions and division order. NumPy-valued
// denominators/custom inputs retain their Python arithmetic instead of being
// silently widened. No peak, match or query result survives this call.
PyObject* quality_overrides(PyObject* module, PyObject* args) {
    PyObject *matches,*precursor,*base;
    double exclusion;
    if (!PyArg_ParseTuple(args,"OOOd",&matches,&precursor,&base,&exclusion)) return nullptr;
    if ((!PyList_CheckExact(matches) && !PyTuple_CheckExact(matches)) ||
            !std::isfinite(exclusion)) Py_RETURN_NONE;
    double precursor_mz, denominator;
    if (builtin_number(precursor,precursor_mz)!=1 || builtin_number(base,denominator)!=1 || denominator<=0) {
        if (PyErr_Occurred()) return nullptr;
        Py_RETURN_NONE;
    }
    auto* state=static_cast<ModuleState*>(PyModule_GetState(module));
    const auto mz_offset=active_offset(state,0), relative_offset=active_offset(state,1);
    const auto fragment_offset=active_offset(state,2), peak_offset=active_offset(state,3);
    const auto intensity_offset=active_offset(state,4);
    Py_ssize_t fragment_offsets[6];
    active_offsets(state,state->fragment_type,state->fragment_fields,state->fragment_descriptors,
                   state->fragment_offsets,fragment_offsets);
    Owned result(PyDict_New());if (!result.value) return nullptr;
    for (Py_ssize_t i=0; i<PySequence_Fast_GET_SIZE(matches); ++i) {
        auto* match=PySequence_Fast_GET_ITEM(matches,i);
        if (reinterpret_cast<PyObject*>(Py_TYPE(match))!=state->match_type) Py_RETURN_NONE;
        Owned fragment(peak_or_match_field(match,state,2,fragment_offset));
        Owned peak(peak_or_match_field(match,state,3,peak_offset));
        if (!fragment.value || !peak.value) return nullptr;
        if (reinterpret_cast<PyObject*>(Py_TYPE(fragment.value))!=state->fragment_type ||
                reinterpret_cast<PyObject*>(Py_TYPE(peak.value))!=state->peak_type) Py_RETURN_NONE;
        Owned type(model_field(fragment.value,state->fragment_fields[2],fragment_offsets[2]));
        if (!type.value) return nullptr;
        if (!PyUnicode_CheckExact(type.value)) Py_RETURN_NONE;
        if (PyUnicode_CompareWithASCIIString(type.value,"Precursor Ion")==0) continue;
        Owned mz(peak_or_match_field(peak.value,state,0,mz_offset));if (!mz.value) return nullptr;
        double mass;
        const int supported_mass=pool_number(mz.value,state,mass);
        if (supported_mass<0) return nullptr;if (!supported_mass) Py_RETURN_NONE;
        if (std::abs(mass-precursor_mz)<=exclusion) continue;
        Owned raw(peak_or_match_field(peak.value,state,4,intensity_offset));
        Owned relative(peak_or_match_field(peak.value,state,1,relative_offset));
        if (!raw.value || !relative.value) return nullptr;
        double raw_intensity, relative_intensity;
        const int a=pool_number(raw.value,state,raw_intensity);
        if (a<0) return nullptr;if (!a) Py_RETURN_NONE;
        const int b=pool_number(relative.value,state,relative_intensity);
        if (b<0) return nullptr;if (!b) Py_RETURN_NONE;
        const double adjusted=std::min(raw_intensity/denominator,1.0);
        if (adjusted>relative_intensity) {
            Owned key(PyLong_FromVoidPtr(match)), value(PyFloat_FromDouble(adjusted));
            if (!key.value || !value.value || PyDict_SetItem(result.value,key.value,value.value)<0) return nullptr;
        }
    }
    return Py_NewRef(result.value);
}

PyObject* quality_overrides_query(PyObject* module,PyObject* args) {
    PyObject *matches,*capsule;
    if (!PyArg_ParseTuple(args,"OO",&matches,&capsule)) return nullptr;
    auto* frame=static_cast<QueryFrame*>(PyCapsule_GetPointer(capsule,QUERY_CAPSULE));
    if (!frame) return nullptr;
    if (!PyList_CheckExact(matches) && !PyTuple_CheckExact(matches)) Py_RETURN_NONE;
    auto* state=static_cast<ModuleState*>(PyModule_GetState(module));
    const auto fragment_offset=active_offset(state,2),peak_offset=active_offset(state,3);
    Py_ssize_t fragment_offsets[6];
    active_offsets(state,state->fragment_type,state->fragment_fields,state->fragment_descriptors,state->fragment_offsets,fragment_offsets);
    Owned result(PyDict_New());if (!result.value) return nullptr;
    if (frame->denominator<=0.) return Py_NewRef(result.value);
    for (Py_ssize_t i=0;i<PySequence_Fast_GET_SIZE(matches);++i) {
        auto* match=PySequence_Fast_GET_ITEM(matches,i);
        if (reinterpret_cast<PyObject*>(Py_TYPE(match))!=state->match_type) Py_RETURN_NONE;
        Owned fragment(peak_or_match_field(match,state,2,fragment_offset));
        Owned peak(peak_or_match_field(match,state,3,peak_offset));
        if (!fragment.value || !peak.value) return nullptr;
        if (reinterpret_cast<PyObject*>(Py_TYPE(fragment.value))!=state->fragment_type) Py_RETURN_NONE;
        Owned type(model_field(fragment.value,state->fragment_fields[2],fragment_offsets[2]));if (!type.value) return nullptr;
        if (!PyUnicode_CheckExact(type.value)) Py_RETURN_NONE;
        if (PyUnicode_CompareWithASCIIString(type.value,"Precursor Ion")==0) continue;
        const auto position=frame->position(peak.value);if (position<0) Py_RETURN_NONE;
        if (std::abs(frame->masses[position]-frame->precursor)<=frame->exclusion) continue;
        const double adjusted=std::min(frame->raw[position]/frame->denominator,1.0);
        if (adjusted>frame->relative[position]) {
            Owned key(PyLong_FromVoidPtr(match)),value(PyFloat_FromDouble(adjusted));
            if (!key.value || !value.value || PyDict_SetItem(result.value,key.value,value.value)<0) return nullptr;
        }
    }
    return Py_NewRef(result.value);
}
double saturation(double intensity, double half) {
    intensity = std::min(std::max(intensity, 0.0), 1.0);
    if (intensity <= 0.0) return 0.0;
    return std::min(((1.0 + half) * intensity) / (intensity + half), 1.0);
}
struct Pool {
    Py_ssize_t total{}, count{};
    double weight{}, half{}, quality{};
    long kind{};
    std::vector<std::pair<double,double>> distinct_peaks;
    std::vector<double> group_quality;
};
constexpr const char* POOL_CAPSULE = "lipidgate.pool_evidence_plan.v1";
struct PoolRoute {std::uintptr_t identity;long mask;Py_ssize_t group;std::uint64_t bit;};
struct CountProgram {
    Py_ssize_t cap;std::uint64_t unassigned;
    std::vector<std::pair<std::uint64_t,Py_ssize_t>> groups;
};
// The same prepared-scoring branches expressed on unboxed counts. No formula,
// threshold or class policy is supplied here: Python prepares those values from
// the revision-checked reference rule. Unsupported values retain Python policy.
struct PreparedDecision {
    bool enabled{},positive{},has_fah{},has_tokens{},ratio_uses_fa{},loss_resolves{};
    bool has_profile{},profile_applies{},profile_species_fallback{},confirm{};
    bool allow_hg{},hg_complete{},single_species{},hg_dominant{},support_scoring{};
    bool lyso{},loss_only{},require_loss{};
    Py_ssize_t hg_total{},loss_total{},hg_required{},negative_hg_required{};
    Py_ssize_t fah_required{},chain_fah_required{},loss_required{},ratio_required{};
    Py_ssize_t profile_required{},expected_chain_count{},signature_required{},fah_min_matches{};
    double fah_min_pool{},fah_hg_dominant_min_pool{},fah_min_quality{},hg_min_relative{};
};
bool read_decision(PyObject* config, PreparedDecision& d) {
    if (!PyDict_CheckExact(config)) return false;
    auto integer=[&](const char* key,Py_ssize_t& value) {
        auto* item=PyDict_GetItemString(config,key);
        if (!item || (!PyLong_CheckExact(item) && !PyBool_Check(item))) return false;
        value=PyLong_AsSsize_t(item);
        if (PyErr_Occurred()) {PyErr_Clear();return false;}
        return value>=0 && value<=1000000;
    };
    auto flag=[&](const char* key,bool& value) {
        Py_ssize_t n;if (!integer(key,n) || n>1) return false;
        value=n!=0;return true;
    };
    auto number=[&](const char* key,double& value) {
        auto* item=PyDict_GetItemString(config,key);
        if (!item || (!PyFloat_CheckExact(item) && !PyLong_CheckExact(item))) return false;
        value=PyFloat_AsDouble(item);
        if (PyErr_Occurred()) {PyErr_Clear();return false;}
        return std::isfinite(value);
    };
#define FLAG(name) if (!flag(#name,d.name)) return false
#define COUNT(name) if (!integer(#name,d.name)) return false
#define NUMBER(name) if (!number(#name,d.name)) return false
    FLAG(positive);FLAG(has_fah);FLAG(has_tokens);FLAG(ratio_uses_fa);FLAG(loss_resolves);
    FLAG(has_profile);FLAG(profile_applies);FLAG(profile_species_fallback);FLAG(confirm);
    FLAG(allow_hg);FLAG(hg_complete);FLAG(single_species);FLAG(hg_dominant);FLAG(support_scoring);
    FLAG(lyso);FLAG(loss_only);FLAG(require_loss);
    COUNT(hg_total);COUNT(loss_total);COUNT(hg_required);COUNT(negative_hg_required);
    COUNT(fah_required);COUNT(chain_fah_required);COUNT(loss_required);COUNT(ratio_required);
    COUNT(profile_required);COUNT(expected_chain_count);COUNT(signature_required);COUNT(fah_min_matches);
    NUMBER(fah_min_pool);NUMBER(fah_hg_dominant_min_pool);NUMBER(fah_min_quality);NUMBER(hg_min_relative);
#undef FLAG
#undef COUNT
#undef NUMBER
    d.enabled=true;return true;
}
enum Resolution {CLASS_LEVEL, SPECIES_LEVEL, CHAIN_LEVEL, TENTATIVE_SPECIES, TENTATIVE_CHAIN};
enum Reason {NO_REASON, CHAIN_CONFIRMATION, SINGLE_CHAIN, MISSING_CHAIN_INFO, MISSING_CHAIN_ANNOTATION,
             LYSO_HG_ONLY, MISSING_POSITIVE_GATE, FAH_FALLBACK, HG_FALLBACK, MISSING_REQUIRED};
enum Missing {MISSING_HG=1,MISSING_LOSS=2,MISSING_FAH=4,MISSING_CONFIRMATION=8};
struct Decision {int missing{},resolution{},reason{};bool reject{};};
Decision decide(const PreparedDecision& p,const std::array<Py_ssize_t,5>& counts,
                Py_ssize_t hg_count,bool has_chain,Py_ssize_t chain_peaks,Py_ssize_t observed,
                Py_ssize_t signature_count,bool raw_loss,double hg_max,const Pool& fah,const Pool& hg) {
    const auto fah_count=counts[0],loss_count=counts[1],fa_loss_count=counts[2];
    const auto chain_loss_count=counts[3],profile_count=counts[4];
    Decision d;
    if (p.positive) {
        if (p.hg_total>0 && hg_count<p.hg_required) d.missing|=MISSING_HG;
        if (!p.hg_total && p.loss_total>0 && loss_count<p.loss_required) d.missing|=MISSING_LOSS;
        if (p.hg_total>0 && p.loss_total>0 && p.ratio_required>0 &&
                (p.ratio_uses_fa ? fa_loss_count : loss_count)<p.ratio_required) d.missing|=MISSING_LOSS;
        if (!p.hg_total && !p.loss_total && p.has_fah && fah_count<p.fah_required) d.missing|=MISSING_FAH;
    } else {
        if (p.has_fah && fah_count<p.fah_required) d.missing|=MISSING_FAH;
        if (p.hg_total>0 && hg_count<p.negative_hg_required) d.missing|=MISSING_HG;
        if ((p.loss_only || p.require_loss) && loss_count<1) d.missing|=MISSING_LOSS;
    }
    if (d.missing && (p.positive ? (d.missing!=MISSING_HG && (d.missing&MISSING_HG)) : d.missing!=MISSING_HG)) {
        d.reject=true;return d;
    }
    int chain_resolution=-1,chain_reason=NO_REASON;
    if (p.confirm && has_chain) {
        if (chain_peaks<2) {chain_resolution=TENTATIVE_SPECIES;chain_reason=CHAIN_CONFIRMATION;}
        else {chain_resolution=CHAIN_LEVEL;chain_reason=(observed<2 || observed<p.expected_chain_count) ? SINGLE_CHAIN : NO_REASON;}
    }
    int configured=chain_resolution,configured_reason=chain_reason;
    if (configured<0 && p.has_profile) {
        if (p.profile_required>0 && profile_count>=p.profile_required) configured=CHAIN_LEVEL;
        else if ((p.profile_required>0 || p.profile_applies) && p.profile_species_fallback) {
            configured=SPECIES_LEVEL;configured_reason=MISSING_CHAIN_INFO;
        }
    }
    if (!p.has_tokens) {d.resolution=SPECIES_LEVEL;d.reason=MISSING_CHAIN_ANNOTATION;}
    else if (configured==CHAIN_LEVEL) {d.resolution=configured;d.reason=configured_reason;}
    else if (p.positive) {
        if ((p.has_fah && fah_count>=p.chain_fah_required) ||
                (p.loss_resolves && chain_loss_count>=p.loss_required)) d.resolution=CHAIN_LEVEL;
        else if (hg_count>=1) {
            if (configured>=0) {d.resolution=configured;d.reason=configured_reason;}
            else if (!p.has_fah && !p.loss_total) {
                if (!p.allow_hg || p.hg_complete || (p.signature_required>0 && signature_count>=p.signature_required))
                    d.resolution=CHAIN_LEVEL;
                else {d.resolution=CLASS_LEVEL;d.reason=LYSO_HG_ONLY;}
            } else {d.resolution=SPECIES_LEVEL;d.reason=MISSING_CHAIN_INFO;}
        } else if (loss_count>=1) {d.resolution=SPECIES_LEVEL;d.reason=MISSING_CHAIN_INFO;}
        else {d.resolution=CLASS_LEVEL;d.reason=MISSING_POSITIVE_GATE;}
    } else if (!p.has_fah && p.loss_only) {
        d.resolution=raw_loss ? CHAIN_LEVEL : SPECIES_LEVEL;d.reason=raw_loss ? NO_REASON : MISSING_CHAIN_INFO;
    } else {
        d.resolution=fah_count>=p.chain_fah_required ? CHAIN_LEVEL : SPECIES_LEVEL;
        d.reason=d.resolution==CHAIN_LEVEL ? NO_REASON : MISSING_CHAIN_INFO;
    }
    const bool fah_low=(p.positive && d.missing==MISSING_HG && fah.total>0 && hg.total>0 && hg.count==0 &&
        fah.count>=p.fah_min_matches && (!p.support_scoring || fah.quality>=p.fah_min_quality) &&
        fah.weight*fah.quality>=(p.hg_dominant ? p.fah_hg_dominant_min_pool : p.fah_min_pool));
    const bool hg_low=(d.resolution!=CHAIN_LEVEL && p.hg_dominant && !(d.missing&MISSING_HG) &&
        hg.total>0 && hg.count>0 && hg_max>=p.hg_min_relative);
    const bool weak_hg=(d.resolution!=CHAIN_LEVEL && !d.missing && p.hg_dominant && hg_count>0 && !hg_low);
    if (fah_low) {d.missing=0;d.resolution=TENTATIVE_CHAIN;d.reason=FAH_FALLBACK;}
    else if (hg_low) {d.missing=0;d.resolution=TENTATIVE_SPECIES;d.reason=HG_FALLBACK;}
    else if (weak_hg) {d.missing=MISSING_HG;d.resolution=CLASS_LEVEL;d.reason=MISSING_REQUIRED;}
    else if (d.missing) {d.resolution=CLASS_LEVEL;d.reason=MISSING_REQUIRED;}
    else if (p.lyso && !p.has_fah) {
        if (d.resolution==CHAIN_LEVEL) d.reason=NO_REASON;
        else if (p.single_species) {d.resolution=SPECIES_LEVEL;d.reason=LYSO_HG_ONLY;}
        else if (p.hg_complete && hg_count>=p.hg_required) {d.resolution=CHAIN_LEVEL;d.reason=NO_REASON;}
        else {d.resolution=CLASS_LEVEL;d.reason=LYSO_HG_ONLY;}
    }
    if (!d.missing && chain_resolution>=0 && chain_reason) {
        d.resolution=chain_resolution;d.reason=chain_reason;
        if (chain_reason==CHAIN_CONFIRMATION) d.missing=MISSING_CONFIRMATION;
    }
    return d;
}
struct PreparedPool {
    std::array<Pool,4> pools;
    Py_ssize_t count{},groups{};int identical{};double replicate{};
    std::vector<PoolRoute> routes;
    std::vector<PoolRoute> occurrences;
    std::uint64_t precursor{};
    bool evidence{};
    std::uint64_t hg{},chain{},raw_loss{},signature{};
    std::vector<std::uint64_t> expected_chains;
    std::vector<CountProgram> programs;
    PreparedDecision decision;
};
void release_pool(PyObject* capsule) {
    auto* value=static_cast<PreparedPool*>(PyCapsule_GetPointer(capsule,POOL_CAPSULE));
    if (value) delete value;else PyErr_Clear();
}
PyObject* create_pool_plan(PyObject*,PyObject* args) {
    PyObject *routes,*configuration,*groups,*statistics=Py_None,*decision=Py_None;
    if (!PyArg_ParseTuple(args,"OOO|OO",&routes,&configuration,&groups,&statistics,&decision)) return nullptr;
    if (!PyDict_CheckExact(routes) || !PyTuple_CheckExact(configuration) ||
            !PyTuple_CheckExact(groups) || (PyTuple_GET_SIZE(groups)!=2 && PyTuple_GET_SIZE(groups)!=3)) Py_RETURN_NONE;
    try {
        auto plan=std::make_unique<PreparedPool>();
        plan->count=PyTuple_GET_SIZE(configuration);
        if (plan->count<1 || plan->count>4) Py_RETURN_NONE;
        plan->groups=PyLong_AsSsize_t(PyTuple_GET_ITEM(groups,0));
        plan->identical=PyObject_IsTrue(PyTuple_GET_ITEM(groups,1));
        plan->replicate=PyTuple_GET_SIZE(groups)==3 ? PyFloat_AsDouble(PyTuple_GET_ITEM(groups,2)) : 0.5;
        if (PyErr_Occurred()) return nullptr;
        if (plan->groups<0 || plan->groups>4096 || !std::isfinite(plan->replicate)) Py_RETURN_NONE;
        for (Py_ssize_t i=0;i<plan->count;++i) {
            auto* item=PyTuple_GET_ITEM(configuration,i);
            if (!PyTuple_CheckExact(item) || PyTuple_GET_SIZE(item)!=5) Py_RETURN_NONE;
            auto& pool=plan->pools[i];
            pool.total=PyLong_AsSsize_t(PyTuple_GET_ITEM(item,1));
            pool.weight=PyFloat_AsDouble(PyTuple_GET_ITEM(item,2));
            pool.half=PyFloat_AsDouble(PyTuple_GET_ITEM(item,3));
            pool.kind=PyLong_AsLong(PyTuple_GET_ITEM(item,4));
            if (PyErr_Occurred()) return nullptr;
            if (pool.total<0 || pool.kind<0 || pool.kind>3 || pool.half<0 ||
                    !std::isfinite(pool.half) || !std::isfinite(pool.weight)) Py_RETURN_NONE;
        }
        PyObject* bits=Py_None;
        if (statistics!=Py_None) {
            if (!PyTuple_CheckExact(statistics) || (PyTuple_GET_SIZE(statistics)!=7 && PyTuple_GET_SIZE(statistics)!=8)) Py_RETURN_NONE;
            auto* programs=PyTuple_GET_ITEM(statistics,0);bits=PyTuple_GET_ITEM(statistics,1);
            auto* expected=PyTuple_GET_ITEM(statistics,4);
            if (!PyTuple_CheckExact(programs) || !PyDict_CheckExact(bits) || !PyTuple_CheckExact(expected)) Py_RETURN_NONE;
            plan->hg=PyLong_AsUnsignedLongLong(PyTuple_GET_ITEM(statistics,2));
            plan->chain=PyLong_AsUnsignedLongLong(PyTuple_GET_ITEM(statistics,3));
            plan->raw_loss=PyLong_AsUnsignedLongLong(PyTuple_GET_ITEM(statistics,5));
            plan->signature=PyLong_AsUnsignedLongLong(PyTuple_GET_ITEM(statistics,6));
            if (PyTuple_GET_SIZE(statistics)==8)
                plan->precursor=PyLong_AsUnsignedLongLong(PyTuple_GET_ITEM(statistics,7));
            if (PyErr_Occurred()) return nullptr;
            for (Py_ssize_t i=0;i<PyTuple_GET_SIZE(expected);++i) {
                plan->expected_chains.push_back(PyLong_AsUnsignedLongLong(PyTuple_GET_ITEM(expected,i)));
                if (PyErr_Occurred()) return nullptr;
            }
            for (Py_ssize_t i=0;i<PyTuple_GET_SIZE(programs);++i) {
                auto* item=PyTuple_GET_ITEM(programs,i);
                if (!PyTuple_CheckExact(item) || PyTuple_GET_SIZE(item)!=3) Py_RETURN_NONE;
                CountProgram program;
                program.cap=PyLong_AsSsize_t(PyTuple_GET_ITEM(item,0));
                program.unassigned=PyLong_AsUnsignedLongLong(PyTuple_GET_ITEM(item,1));
                auto* pairs=PyTuple_GET_ITEM(item,2);
                if (PyErr_Occurred()) return nullptr;
                if (program.cap< -1 || !PyTuple_CheckExact(pairs)) Py_RETURN_NONE;
                for (Py_ssize_t j=0;j<PyTuple_GET_SIZE(pairs);++j) {
                    auto* pair=PyTuple_GET_ITEM(pairs,j);
                    if (!PyTuple_CheckExact(pair) || PyTuple_GET_SIZE(pair)!=2) Py_RETURN_NONE;
                    auto mask=PyLong_AsUnsignedLongLong(PyTuple_GET_ITEM(pair,0));
                    auto weight=PyLong_AsSsize_t(PyTuple_GET_ITEM(pair,1));
                    if (PyErr_Occurred()) return nullptr;
                    if (weight<0 || weight>1000000) Py_RETURN_NONE;
                    program.groups.emplace_back(mask,weight);
                }
                plan->programs.push_back(std::move(program));
            }
            plan->evidence=true;
        }
        if (decision!=Py_None && plan->programs.size()==5 && plan->count>=2 &&
                PyUnicode_CheckExact(PyTuple_GET_ITEM(PyTuple_GET_ITEM(configuration,0),0)) &&
                PyUnicode_CheckExact(PyTuple_GET_ITEM(PyTuple_GET_ITEM(configuration,1),0)) &&
                PyUnicode_CompareWithASCIIString(PyTuple_GET_ITEM(PyTuple_GET_ITEM(configuration,0),0),"fah")==0 &&
                PyUnicode_CompareWithASCIIString(PyTuple_GET_ITEM(PyTuple_GET_ITEM(configuration,1),0),"hg")==0)
            read_decision(decision,plan->decision);
        Py_ssize_t pos=0;PyObject *key,*value;
        plan->routes.reserve(PyDict_Size(routes));
        while (PyDict_Next(routes,&pos,&key,&value)) {
            if (!PyLong_CheckExact(key) || !PyTuple_CheckExact(value) || PyTuple_GET_SIZE(value)!=2) Py_RETURN_NONE;
            auto* identity=PyLong_AsVoidPtr(key);
            const auto mask=PyLong_AsLong(PyTuple_GET_ITEM(value,0));
            const auto group=PyLong_AsSsize_t(PyTuple_GET_ITEM(value,1));
            if (PyErr_Occurred()) return nullptr;
            if (mask<0 || mask>=(1L<<plan->count) || group< -1 || group>=plan->groups) Py_RETURN_NONE;
            std::uint64_t bit=0;
            if (bits!=Py_None) {
                auto* found=PyDict_GetItemWithError(bits,key);
                if (!found) {if (PyErr_Occurred()) return nullptr;Py_RETURN_NONE;}
                bit=PyLong_AsUnsignedLongLong(found);if (PyErr_Occurred()) return nullptr;
            }
            plan->routes.push_back({reinterpret_cast<std::uintptr_t>(identity),mask,group,bit});
        }
        std::sort(plan->routes.begin(),plan->routes.end(),[](const auto& a,const auto& b){return a.identity<b.identity;});
        if (statistics!=Py_None && PyTuple_GET_SIZE(statistics)==8 && plan->routes.size()<=64) {
            plan->occurrences.resize(plan->routes.size());
            std::uint64_t seen=0;bool valid=true;
            for (const auto& route:plan->routes) {
                if (!route.bit || (route.bit & (route.bit-1)) || (seen & route.bit)) {valid=false;break;}
                const auto index=__builtin_ctzll(route.bit);
                if (index>=plan->occurrences.size()) {valid=false;break;}
                seen |= route.bit;plan->occurrences[index]=route;
            }
            if (!valid) plan->occurrences.clear();
        }
        auto* capsule=PyCapsule_New(plan.get(),POOL_CAPSULE,release_pool);
        if (capsule) plan.release();return capsule;
    } catch (const std::bad_alloc&) {return PyErr_NoMemory();}
      catch (...) {PyErr_SetString(PyExc_RuntimeError,"Native pool preparation failed");return nullptr;}
}
PyObject* pool_metrics_impl(PyObject* module, PyObject* args, bool execute_policy, bool normalize=false, bool indexed=false) {
    PyObject *matches, *routes=nullptr, *configuration=nullptr, *groups=nullptr, *overrides=Py_None, *normalization=nullptr;
    PyObject* query_object=nullptr;QueryFrame* frame=nullptr;
    PreparedPool* prepared=nullptr;
    if (PyTuple_GET_SIZE(args)==4 && !normalize) {
        PyObject* capsule;
        if (!PyArg_ParseTuple(args,"OOOO",&matches,&capsule,&overrides,&query_object)) return nullptr;
        prepared=static_cast<PreparedPool*>(PyCapsule_GetPointer(capsule,POOL_CAPSULE));
        if (!prepared) return nullptr;
    } else if (PyTuple_GET_SIZE(args)==3) {
        PyObject* capsule;
        if (!PyArg_ParseTuple(args,"OOO",&matches,&capsule,normalize ? &normalization : &overrides)) return nullptr;
        prepared=static_cast<PreparedPool*>(PyCapsule_GetPointer(capsule,POOL_CAPSULE));
        if (!prepared) return nullptr;
    } else if (!PyArg_ParseTuple(args,"OOOOO",&matches,&routes,&configuration,&groups,&overrides)) return nullptr;
    if (execute_policy && (!prepared || !prepared->decision.enabled)) Py_RETURN_NONE;
    if ((!PyList_CheckExact(matches) && !PyTuple_CheckExact(matches)) || (!prepared && (!PyDict_CheckExact(routes) ||
            !PyTuple_CheckExact(configuration) || !PyTuple_CheckExact(groups) ||
            (PyTuple_GET_SIZE(groups)!=2 && PyTuple_GET_SIZE(groups)!=3))) ||
            (overrides!=Py_None && !PyDict_CheckExact(overrides))) Py_RETURN_NONE;
    const auto count = PySequence_Fast_GET_SIZE(matches), pool_count = prepared ? prepared->count : PyTuple_GET_SIZE(configuration);
    if (pool_count<1 || pool_count>4) Py_RETURN_NONE;
    double precursor_mz=0., denominator=0., exclusion=0.;
    if (normalize && normalization && PyCapsule_IsValid(normalization,QUERY_CAPSULE)) query_object=normalization;
    if (query_object) {
        frame=static_cast<QueryFrame*>(PyCapsule_GetPointer(query_object,QUERY_CAPSULE));
        if (!frame) return nullptr;
        precursor_mz=frame->precursor;denominator=frame->denominator;exclusion=frame->exclusion;
    }
    if (indexed && (!frame || !prepared || count!=static_cast<Py_ssize_t>(prepared->occurrences.size()) || !prepared->evidence)) Py_RETURN_NONE;
    if (normalize) {
        if (!frame) {
            if (!prepared || !normalization || !PyTuple_CheckExact(normalization) || PyTuple_GET_SIZE(normalization)!=3) Py_RETURN_NONE;
            const int a=builtin_number(PyTuple_GET_ITEM(normalization,0),precursor_mz);
            const int b=builtin_number(PyTuple_GET_ITEM(normalization,1),denominator);
            const int c=builtin_number(PyTuple_GET_ITEM(normalization,2),exclusion);
            if (a<0 || b<0 || c<0) return nullptr;
            if (!a || !b || !c || denominator<=0.) Py_RETURN_NONE;
        }
    }
    auto* state = static_cast<ModuleState*>(PyModule_GetState(module));
    const auto mass_offset=indexed ? -1 : active_offset(state,0), intensity_offset=indexed ? -1 : active_offset(state,1);
    const auto fragment_offset=indexed ? -1 : active_offset(state,2), peak_offset=indexed ? -1 : active_offset(state,3);
    const auto raw_intensity_offset=normalize && !indexed ? active_offset(state,4) : -1;
    Py_ssize_t fragment_offsets[6];
    if (normalize && !indexed) active_offsets(state,state->fragment_type,state->fragment_fields,state->fragment_descriptors,
                                  state->fragment_offsets,fragment_offsets);
    const auto group_count = prepared ? prepared->groups : PyLong_AsSsize_t(PyTuple_GET_ITEM(groups,0));
    const int identical = prepared ? prepared->identical : PyObject_IsTrue(PyTuple_GET_ITEM(groups,1));
    const double replicate = prepared ? prepared->replicate : PyTuple_GET_SIZE(groups)==3 ? PyFloat_AsDouble(PyTuple_GET_ITEM(groups,2)) : 0.5;
    if (PyErr_Occurred()) return nullptr;
    if (group_count<0 || group_count>4096 || !std::isfinite(replicate)) Py_RETURN_NONE;
    try {
        std::array<Pool,4> pools;
        for (Py_ssize_t i=0; i<pool_count; ++i) {
            if (prepared) {
                pools[i]=prepared->pools[i];
                if (pools[i].kind==3) pools[i].group_quality.resize(group_count,0.0);
                continue;
            }
            auto* item = PyTuple_GET_ITEM(configuration,i);
            if (!PyTuple_CheckExact(item) || PyTuple_GET_SIZE(item)!=5) Py_RETURN_NONE;
            auto& pool = pools[i];
            pool.total = PyLong_AsSsize_t(PyTuple_GET_ITEM(item,1));
            pool.weight = PyFloat_AsDouble(PyTuple_GET_ITEM(item,2));
            pool.half = PyFloat_AsDouble(PyTuple_GET_ITEM(item,3));
            pool.kind = PyLong_AsLong(PyTuple_GET_ITEM(item,4));
            if (PyErr_Occurred()) return nullptr;
            if (pool.total<0 || !std::isfinite(pool.weight) || !std::isfinite(pool.half) ||
                    pool.half<0 || pool.kind<0 || pool.kind>3) Py_RETURN_NONE;
            if (pool.kind==3) pool.group_quality.resize(group_count,0.0);
        }
        std::uint64_t observed_mask=0;double hg_max=0.;bool has_hg=false;
        std::vector<double> chain_peaks;
        const bool has_overrides=overrides!=Py_None && PyDict_Size(overrides)>0;
        auto** items = PySequence_Fast_ITEMS(matches);
        for (Py_ssize_t j=0; j<count; ++j) {
            auto* match = items[j];
            if (!indexed && reinterpret_cast<PyObject*>(Py_TYPE(match))!=state->match_type) Py_RETURN_NONE;
            Owned fragment(indexed ? nullptr : peak_or_match_field(match,state,2,fragment_offset));
            Owned peak(indexed ? nullptr : peak_or_match_field(match,state,3,peak_offset));
            if (!indexed && (!fragment.value || !peak.value)) return nullptr;
            if (!indexed && (reinterpret_cast<PyObject*>(Py_TYPE(fragment.value))!=state->fragment_type ||
                    reinterpret_cast<PyObject*>(Py_TYPE(peak.value))!=state->peak_type)) Py_RETURN_NONE;
            Py_ssize_t peak_position=-1;
            if (indexed) {
                if (!PyLong_CheckExact(match)) Py_RETURN_NONE;
                peak_position=PyLong_AsSsize_t(match);
                if (PyErr_Occurred()) {PyErr_Clear();Py_RETURN_NONE;}
                if (peak_position<0) continue;
                if (peak_position>=static_cast<Py_ssize_t>(frame->masses.size())) Py_RETURN_NONE;
            } else if (frame) peak_position=frame->position(peak.value);
            if (frame && peak_position<0) Py_RETURN_NONE;
            long mask;Py_ssize_t group;std::uint64_t evidence_bit=0;
            if (indexed) {
                const auto& route=prepared->occurrences[j];
                mask=route.mask;group=route.group;evidence_bit=route.bit;observed_mask|=evidence_bit;
            } else if (prepared) {
                const auto identity=reinterpret_cast<std::uintptr_t>(fragment.value);
                auto found=std::lower_bound(prepared->routes.begin(),prepared->routes.end(),identity,
                    [](const auto& route,auto key){return route.identity<key;});
                if (found==prepared->routes.end() || found->identity!=identity) Py_RETURN_NONE;
                mask=found->mask;group=found->group;evidence_bit=found->bit;observed_mask|=evidence_bit;
            } else {
                Owned key(PyLong_FromVoidPtr(fragment.value)); if (!key.value) return nullptr;
                auto* route = PyDict_GetItemWithError(routes,key.value);
                if (!route) { if (PyErr_Occurred()) return nullptr; Py_RETURN_NONE; }
                if (!PyTuple_CheckExact(route) || PyTuple_GET_SIZE(route)!=2) Py_RETURN_NONE;
                mask=PyLong_AsLong(PyTuple_GET_ITEM(route,0));group=PyLong_AsSsize_t(PyTuple_GET_ITEM(route,1));
                if (PyErr_Occurred()) return nullptr;
                if (mask<0 || mask>=(1L<<pool_count) || group< -1 || group>=group_count) Py_RETURN_NONE;
            }
            double relative=0.0;
            Owned match_id(has_overrides ? PyLong_FromVoidPtr(match) : nullptr);
            if (has_overrides && !match_id.value) return nullptr;
            auto* override = has_overrides ? PyDict_GetItemWithError(overrides,match_id.value) : nullptr;
            if (!override && PyErr_Occurred()) return nullptr;
            Owned original(override || frame ? nullptr : peak_or_match_field(peak.value,state,1,intensity_offset));
            if (!override && !frame && !original.value) return nullptr;
            const int supported = frame && !override ? (relative=frame->relative[peak_position],1) :
                pool_number(override ? override : original.value,state,relative);
            if (supported<0) return nullptr;
            if (!supported) Py_RETURN_NONE;
            const double raw_relative=frame ? frame->relative[peak_position] : relative;
            if (normalize && denominator>0.) {
                Owned type(indexed ? nullptr : model_field(fragment.value,state->fragment_fields[2],fragment_offsets[2]));
                if (!indexed && !type.value) return nullptr;
                if (!indexed && !PyUnicode_CheckExact(type.value)) Py_RETURN_NONE;
                if (indexed ? !(evidence_bit & prepared->precursor) : PyUnicode_CompareWithASCIIString(type.value,"Precursor Ion")!=0) {
                    Owned mz(frame ? nullptr : peak_or_match_field(peak.value,state,0,mass_offset));
                    if (!frame && !mz.value) return nullptr;
                    double mass;const int a=frame ? (mass=frame->masses[peak_position],1) : pool_number(mz.value,state,mass);
                    if (a<0) return nullptr;if (!a) Py_RETURN_NONE;
                    if (std::abs(mass-precursor_mz)>exclusion) {
                        Owned raw(frame ? nullptr : peak_or_match_field(peak.value,state,4,raw_intensity_offset));
                        if (!frame && !raw.value) return nullptr;
                        double intensity;const int b=frame ? (intensity=frame->raw[peak_position],1) : pool_number(raw.value,state,intensity);
                        if (b<0) return nullptr;if (!b) Py_RETURN_NONE;
                        const double adjusted=std::min(intensity/denominator,1.0);
                        if (adjusted>relative) relative=adjusted;
                    }
                }
            }
            if (prepared && prepared->evidence) {
                if (evidence_bit & prepared->hg) {
                    double hg_relative=raw_relative;
                    if (override && !frame) {
                        Owned raw(peak_or_match_field(peak.value,state,1,intensity_offset));if (!raw.value) return nullptr;
                        const int valid=pool_number(raw.value,state,hg_relative);
                        if (valid<0) return nullptr;if (!valid) Py_RETURN_NONE;
                    }
                    if (!has_hg || hg_relative>hg_max) hg_max=hg_relative;
                    has_hg=true;
                }
                if (evidence_bit & prepared->chain) {
                    Owned mz(frame ? nullptr : peak_or_match_field(peak.value,state,0,mass_offset));
                    if (!frame && !mz.value) return nullptr;
                    double mass;const int valid=frame ? (mass=frame->masses[peak_position],1) : pool_number(mz.value,state,mass);
                    if (valid<0) return nullptr;if (!valid) Py_RETURN_NONE;
                    if (std::find(chain_peaks.begin(),chain_peaks.end(),mass)==chain_peaks.end()) chain_peaks.push_back(mass);
                }
            }
            relative = std::min(std::max(relative,0.0),1.0);
            for (Py_ssize_t i=0; i<pool_count; ++i) {
                if (!(mask & (1L<<i))) continue;
                auto& pool = pools[i]; ++pool.count;
                if (!pool.kind) continue;
                const double quality = saturation(relative,pool.half);
                if (pool.kind==1) pool.quality = std::max(pool.quality,quality);
                else if (pool.kind==3) {
                    if (group>=0) pool.group_quality[group]=std::max(pool.group_quality[group],quality);
                } else {
                    Owned mz(frame ? nullptr : peak_or_match_field(peak.value,state,0,mass_offset));
                    if (!frame && !mz.value) return nullptr;
                    double mass;
                    const int valid=frame ? (mass=frame->masses[peak_position],1) : pool_number(mz.value,state,mass);
                    if (valid<0) return nullptr;
                    if (!valid) Py_RETURN_NONE;
                    auto found=std::find_if(pool.distinct_peaks.begin(),pool.distinct_peaks.end(),
                        [mass](const auto& p){return p.first==mass;});
                    if (found==pool.distinct_peaks.end()) pool.distinct_peaks.emplace_back(mass,quality);
                    else found->second=std::max(found->second,quality);
                }
            }
        }
        if (indexed && !observed_mask) Py_RETURN_FALSE;
        Owned result(PyTuple_New(pool_count)); if (!result.value) return nullptr;
        for (Py_ssize_t i=0; i<pool_count; ++i) {
            auto& pool=pools[i];
            if (!pool.kind) {
                const auto target=std::min<Py_ssize_t>(pool.total,3);
                pool.quality=target ? std::min(static_cast<double>(pool.count)/target,1.0) : 0.0;
            } else if (pool.kind==2) {
                double first=0.0, second=0.0;
                for (const auto& pair:pool.distinct_peaks) {
                    const auto q=pair.second;
                    if (q>first) {second=first;first=q;} else if (q>second) second=q;
                }
                pool.quality=(first+second)/2.0;
            } else if (pool.kind==3 && group_count) {
                auto& values=pool.group_quality;
                std::sort(values.begin(),values.end(),std::greater<double>());
                pool.quality=identical ? (values[0]+values[0]*replicate)/2.0 :
                    (values[0]+(group_count>1 ? values[1] : 0.0))/std::min<Py_ssize_t>(2,group_count);
            }
            auto* row=Py_BuildValue("(nddd)",pool.count,
                pool.total ? static_cast<double>(pool.count)/pool.total : 0.0,
                pool.quality,pool.weight*pool.quality);
            if (!row) return nullptr;
            PyTuple_SET_ITEM(result.value,i,row);
        }
        if (execute_policy) {
            std::array<Py_ssize_t,5> counts{};
            for (std::size_t i=0;i<counts.size();++i) {
                const auto& program=prepared->programs[i];
                auto& value=counts[i];value=__builtin_popcountll(observed_mask & program.unassigned);
                for (const auto& group:program.groups) if (observed_mask & group.first) value+=group.second;
                if (program.cap>=0) value=std::min(value,program.cap);
            }
            Py_ssize_t expected=0;for (auto mask:prepared->expected_chains) expected+=!!(observed_mask & mask);
            const auto decision=decide(prepared->decision,counts,
                __builtin_popcountll(observed_mask & prepared->hg),!!(observed_mask & prepared->chain),
                chain_peaks.size(),expected,__builtin_popcountll(observed_mask & prepared->signature),
                !!(observed_mask & prepared->raw_loss),hg_max,pools[0],pools[1]);
            if (decision.reject) Py_RETURN_FALSE;
            return Py_BuildValue("(Oiii)",result.value,decision.missing,decision.resolution,decision.reason);
        }
        if (prepared && prepared->evidence) {
            Owned counts(PyTuple_New(prepared->programs.size()));if (!counts.value) return nullptr;
            for (std::size_t i=0;i<prepared->programs.size();++i) {
                const auto& program=prepared->programs[i];
                Py_ssize_t value=__builtin_popcountll(observed_mask & program.unassigned);
                for (const auto& group:program.groups) if (observed_mask & group.first) value+=group.second;
                if (program.cap>=0) value=std::min(value,program.cap);
                auto* item=PyLong_FromSsize_t(value);if (!item) return nullptr;PyTuple_SET_ITEM(counts.value,i,item);
            }
            Py_ssize_t expected=0;for (auto mask:prepared->expected_chains) expected+=!!(observed_mask & mask);
            Owned evidence(Py_BuildValue("(KOnnndnn)",static_cast<unsigned long long>(observed_mask),counts.value,
                static_cast<Py_ssize_t>(__builtin_popcountll(observed_mask & prepared->hg)),
                static_cast<Py_ssize_t>(__builtin_popcountll(observed_mask & prepared->raw_loss)),
                static_cast<Py_ssize_t>(__builtin_popcountll(observed_mask & prepared->signature)),
                hg_max,static_cast<Py_ssize_t>(chain_peaks.size()),expected));
            if (!evidence.value) return nullptr;
            return PyTuple_Pack(2,result.value,evidence.value);
        }
        Py_INCREF(result.value);return result.value;
    } catch (const std::bad_alloc&) { return PyErr_NoMemory(); }
      catch (...) { PyErr_SetString(PyExc_RuntimeError,"Native pool calculation failed");return nullptr; }
}
PyObject* pool_metrics(PyObject* module, PyObject* args) {return pool_metrics_impl(module,args,false);}
PyObject* score_pool_policy(PyObject* module, PyObject* args) {return pool_metrics_impl(module,args,true);}
PyObject* pool_metrics_normalized(PyObject* module, PyObject* args) {return pool_metrics_impl(module,args,false,true);}
PyObject* score_pool_policy_normalized(PyObject* module, PyObject* args) {return pool_metrics_impl(module,args,true,true);}

PyObject* indexed_scoring_compatible(PyObject* module,PyObject* constructor) {
    auto* state=static_cast<ModuleState*>(PyModule_GetState(module));
    if (constructor!=state->match_type || !state->match_initializer ||
            !PyTuple_CheckExact(state->match_initializer) || PyTuple_GET_SIZE(state->match_initializer)!=2) Py_RETURN_FALSE;
    auto* type=reinterpret_cast<PyTypeObject*>(constructor);
    Owned initializer(PyObject_GetAttrString(constructor,"__init__"));if (!initializer.value) return nullptr;
    if (initializer.value!=PyTuple_GET_ITEM(state->match_initializer,0) || !PyFunction_Check(initializer.value) ||
            PyFunction_GET_CODE(initializer.value)!=PyTuple_GET_ITEM(state->match_initializer,1) ||
            type->tp_new!=PyBaseObject_Type.tp_new || type->tp_alloc!=PyType_GenericAlloc || type->tp_finalize ||
            active_offset(state,2)<0 || active_offset(state,3)<0 || active_offset(state,5)<0) Py_RETURN_FALSE;
    Py_ssize_t fragment_offsets[6];
    active_offsets(state,state->fragment_type,state->fragment_fields,state->fragment_descriptors,state->fragment_offsets,fragment_offsets);
    for (auto offset:fragment_offsets) if (offset<0) Py_RETURN_FALSE;
    Py_RETURN_TRUE;
}

PyObject* indexed_pool_result(PyObject* module,PyObject* args,bool policy) {
    Owned result(pool_metrics_impl(module,args,policy,true,true));
    if (!result.value) return nullptr;
    if (result.value==Py_None || result.value==Py_False) return Py_NewRef(result.value);
    auto* indexes=PyTuple_GET_ITEM(args,0);
    auto* frame=static_cast<QueryFrame*>(PyCapsule_GetPointer(PyTuple_GET_ITEM(args,2),QUERY_CAPSULE));
    if (!frame) return nullptr;
    Py_ssize_t count=0;
    for (Py_ssize_t i=0;i<PySequence_Fast_GET_SIZE(indexes);++i)
        count+=PyLong_AsSsize_t(PySequence_Fast_GET_ITEM(indexes,i))>=0;
    if (PyErr_Occurred()) return nullptr;
    Owned raw(PyTuple_New(count)),relative(PyTuple_New(count));
    if (!raw.value || !relative.value) return nullptr;
    Py_ssize_t target=0;
    for (Py_ssize_t i=0;i<PySequence_Fast_GET_SIZE(indexes);++i) {
        const auto position=PyLong_AsSsize_t(PySequence_Fast_GET_ITEM(indexes,i));
        if (position<0) continue;
        PyTuple_SET_ITEM(raw.value,target,Py_NewRef(frame->raw_objects[position]));
        PyTuple_SET_ITEM(relative.value,target++,Py_NewRef(frame->relative_objects[position]));
    }
    auto* sum=PyDict_GetItemString(PyEval_GetBuiltins(),"sum");
    if (!sum || !PyCFunction_Check(sum)) Py_RETURN_NONE;
    Owned a(PyObject_CallOneArg(sum,raw.value)),b(PyObject_CallOneArg(sum,relative.value));
    if (!a.value || !b.value) return nullptr;
    Owned sums(PyTuple_Pack(2,a.value,b.value));if (!sums.value) return nullptr;
    return PyTuple_Pack(2,result.value,sums.value);
}
PyObject* score_pool_indexes(PyObject* module,PyObject* args) {return indexed_pool_result(module,args,true);}
PyObject* pool_metrics_indexes(PyObject* module,PyObject* args) {return indexed_pool_result(module,args,false);}

// One invocation owns the entire shared-chain greedy rerank. Library pool
// programs are borrowed; all experimental state is local and released here.
struct SharedPoolMatch {
    long mask{};Py_ssize_t group{};std::uintptr_t peak{},fragment{};
    double mass{},relative{};bool chain{};unsigned int penalties{};
};
struct SharedPoolCandidate {
    PreparedPool* plan{};std::vector<SharedPoolMatch> matches;
    Py_ssize_t match_count{};double score{},ppm{};bool bonus{};
    PyObject *score_object{},*metrics{},*configuration{},*capsule{},*prepare_arguments{};
    ~SharedPoolCandidate() {
        Py_XDECREF(score_object);Py_XDECREF(metrics);Py_XDECREF(configuration);Py_XDECREF(capsule);
    }
};
bool shared_score(SharedPoolCandidate& candidate,double factor,double score_min,double score_max,
                  PyObject* sum_function,PyObject* round_function) {
    const auto* plan=candidate.plan;
    auto pools=plan->pools;
    for (Py_ssize_t i=0;i<plan->count;++i)
        if (pools[i].kind==3) pools[i].group_quality.assign(plan->groups,0.0);
    for (const auto& match:candidate.matches) {
        double relative=match.relative;
        if (match.penalties) {
            // CPython pow keeps the exact reference arithmetic, including
            // non-default factors. The unpenalized intensity selects winners.
            Owned base(PyFloat_FromDouble(factor)),exponent(PyLong_FromUnsignedLong(match.penalties));
            if (!base.value || !exponent.value) return false;
            Owned power(PyNumber_Power(base.value,exponent.value,Py_None));if (!power.value) return false;
            relative*=PyFloat_AsDouble(power.value);
            if (PyErr_Occurred()) return false;
        }
        relative=std::min(std::max(relative,0.0),1.0);
        for (Py_ssize_t i=0;i<plan->count;++i) {
            if (!(match.mask & (1L<<i))) continue;
            auto& pool=pools[i];++pool.count;
            if (!pool.kind) continue;
            const double quality=saturation(relative,pool.half);
            if (pool.kind==1) pool.quality=std::max(pool.quality,quality);
            else if (pool.kind==3) {
                if (match.group>=0) pool.group_quality[match.group]=std::max(pool.group_quality[match.group],quality);
            } else {
                auto at=std::find_if(pool.distinct_peaks.begin(),pool.distinct_peaks.end(),
                    [&](const auto& pair){return pair.first==match.mass;});
                if (at==pool.distinct_peaks.end()) pool.distinct_peaks.emplace_back(match.mass,quality);
                else at->second=std::max(at->second,quality);
            }
        }
    }
    Owned metrics(PyTuple_New(plan->count)),parts(PyTuple_New(plan->count));
    if (!metrics.value || !parts.value) return false;
    for (Py_ssize_t i=0;i<plan->count;++i) {
        auto& pool=pools[i];
        if (!pool.kind) {
            const auto target=std::min<Py_ssize_t>(pool.total,3);
            pool.quality=target ? std::min(static_cast<double>(pool.count)/target,1.0) : 0.0;
        } else if (pool.kind==2) {
            double first=0.0,second=0.0;
            for (const auto& pair:pool.distinct_peaks) {
                if (pair.second>first) {second=first;first=pair.second;}
                else if (pair.second>second) second=pair.second;
            }
            pool.quality=(first+second)/2.0;
        } else if (pool.kind==3 && plan->groups) {
            auto& values=pool.group_quality;
            std::sort(values.begin(),values.end(),std::greater<double>());
            pool.quality=plan->identical ? (values[0]+values[0]*plan->replicate)/2.0 :
                (values[0]+(plan->groups>1 ? values[1] : 0.0))/std::min<Py_ssize_t>(2,plan->groups);
        }
        auto* row=Py_BuildValue("(nddd)",pool.count,
            pool.total ? static_cast<double>(pool.count)/pool.total : 0.0,
            pool.quality,pool.weight*pool.quality);
        if (!row) return false;
        PyTuple_SET_ITEM(metrics.value,i,row);
        PyTuple_SET_ITEM(parts.value,i,Py_NewRef(PyTuple_GET_ITEM(row,3)));
    }
    // Python 3.12+ uses a compensated float sum. Do not replace it with a
    // C++ accumulation, or replace Python's decimal rounding with std::round.
    Owned total(PyObject_CallOneArg(sum_function,parts.value));if (!total.value) return false;
    double value=PyFloat_AsDouble(total.value);if (PyErr_Occurred()) return false;
    value=std::min(score_max,std::max(score_min,value));
    if (candidate.bonus) value=std::min(score_max,value+10.0);
    Owned number(PyFloat_FromDouble(value)),digits(PyLong_FromLong(4));
    if (!number.value || !digits.value) return false;
    auto* rounded=PyObject_CallFunctionObjArgs(round_function,number.value,digits.value,nullptr);
    if (!rounded) return false;
    candidate.score=PyFloat_AsDouble(rounded);
    if (PyErr_Occurred()) {Py_DECREF(rounded);return false;}
    Py_XSETREF(candidate.score_object,rounded);
    Py_XSETREF(candidate.metrics,Py_NewRef(metrics.value));
    return true;
}
PyObject* rerank_pool_candidates(PyObject* module,PyObject* args) {
    PyObject *specifications,*query,*selected_groups,*prepare_pool;
    double factor,minimum,ppm_limit,score_min,score_max;
    if (!PyArg_ParseTuple(args,"OOOOddddd",&specifications,&query,&selected_groups,&prepare_pool,
                         &factor,&minimum,&ppm_limit,&score_min,&score_max)) return nullptr;
    if (!PyTuple_CheckExact(specifications) || !PyTuple_CheckExact(selected_groups) || !PyCallable_Check(prepare_pool) ||
            !std::isfinite(factor) || factor<0.0 || factor>1.0 ||
            !std::isfinite(minimum) || !std::isfinite(ppm_limit) || ppm_limit<=0.0 ||
            !std::isfinite(score_min) || !std::isfinite(score_max)) Py_RETURN_NONE;
    auto* frame=static_cast<QueryFrame*>(PyCapsule_GetPointer(query,QUERY_CAPSULE));
    if (!frame) return nullptr;
    auto* state=static_cast<ModuleState*>(PyModule_GetState(module));
    const auto fragment_offset=active_offset(state,2),peak_offset=active_offset(state,3);
    if (active_offset(state,0)<0 || active_offset(state,1)<0 || active_offset(state,4)<0) Py_RETURN_NONE;
    Py_ssize_t fragment_offsets[6];
    active_offsets(state,state->fragment_type,state->fragment_fields,state->fragment_descriptors,
                   state->fragment_offsets,fragment_offsets);
    if (fragment_offset<0 || peak_offset<0) Py_RETURN_NONE;
    for (auto offset:fragment_offsets) if (offset<0) Py_RETURN_NONE;
    auto* builtins=PyEval_GetBuiltins();
    auto* sum_function=PyDict_GetItemString(builtins,"sum");
    auto* round_function=PyDict_GetItemString(builtins,"round");
    if (!sum_function || !round_function || !PyCFunction_Check(sum_function) || !PyCFunction_Check(round_function)) Py_RETURN_NONE;
    try {
        std::vector<std::vector<std::uintptr_t>> groups;
        for (Py_ssize_t i=0;i<PyTuple_GET_SIZE(selected_groups);++i) {
            auto* row=PyTuple_GET_ITEM(selected_groups,i);if (!PyTuple_CheckExact(row)) Py_RETURN_NONE;
            std::vector<std::uintptr_t> group;
            for (Py_ssize_t j=0;j<PyTuple_GET_SIZE(row);++j) {
                auto* item=PyTuple_GET_ITEM(row,j);if (!PyLong_CheckExact(item)) Py_RETURN_NONE;
                auto value=reinterpret_cast<std::uintptr_t>(PyLong_AsVoidPtr(item));
                if (PyErr_Occurred()) return nullptr;
                group.push_back(value);
            }
            if (!group.empty()) groups.push_back(std::move(group));
        }
        const auto count=PyTuple_GET_SIZE(specifications);
        std::vector<std::unique_ptr<SharedPoolCandidate>> candidates;
        candidates.reserve(count);
        std::vector<Py_ssize_t> remaining,selected;
        for (Py_ssize_t i=0;i<count;++i) {
            auto* row=PyTuple_GET_ITEM(specifications,i);
            if (!PyTuple_CheckExact(row) || PyTuple_GET_SIZE(row)!=7) Py_RETURN_NONE;
            auto* matches=PyTuple_GET_ITEM(row,0);
            auto* prepare_arguments=PyTuple_GET_ITEM(row,1);
            auto* normalize=PyTuple_GET_ITEM(row,2);
            auto* chains=PyTuple_GET_ITEM(row,3);
            auto* initial=PyTuple_GET_ITEM(row,4);
            auto* ppm=PyTuple_GET_ITEM(row,5);
            auto* bonus=PyTuple_GET_ITEM(row,6);
            if ((!PyList_CheckExact(matches) && !PyTuple_CheckExact(matches)) || !PyTuple_CheckExact(chains) ||
                    !PyBool_Check(normalize) || !PyBool_Check(bonus) || !PyFloat_CheckExact(initial) ||
                    !PyFloat_CheckExact(ppm)) Py_RETURN_NONE;
            auto candidate=std::make_unique<SharedPoolCandidate>();
            candidate->score=PyFloat_AS_DOUBLE(initial);candidate->ppm=PyFloat_AS_DOUBLE(ppm);
            if (!std::isfinite(candidate->score) || !std::isfinite(candidate->ppm)) Py_RETURN_NONE;
            candidate->score_object=Py_NewRef(initial);
            candidate->match_count=PySequence_Fast_GET_SIZE(matches);
            candidate->bonus=bonus==Py_True;
            if (PyTuple_GET_SIZE(chains)) {
                if (!PyTuple_CheckExact(prepare_arguments) || PyTuple_GET_SIZE(prepare_arguments)!=3) Py_RETURN_NONE;
                candidate->prepare_arguments=prepare_arguments;
                std::vector<bool> chain_flags(candidate->match_count,false);
                for (Py_ssize_t j=0;j<PyTuple_GET_SIZE(chains);++j) {
                    auto* item=PyTuple_GET_ITEM(chains,j);if (!PyLong_CheckExact(item)) Py_RETURN_NONE;
                    auto index=PyLong_AsSsize_t(item);if (PyErr_Occurred()) return nullptr;
                    if (index<0 || index>=candidate->match_count) Py_RETURN_NONE;
                    chain_flags[index]=true;
                }
                std::vector<PyObject*> match_identities;
                for (Py_ssize_t j=0;j<candidate->match_count;++j) {
                    auto* match=PySequence_Fast_GET_ITEM(matches,j);
                    if (reinterpret_cast<PyObject*>(Py_TYPE(match))!=state->match_type) Py_RETURN_NONE;
                    if (std::find(match_identities.begin(),match_identities.end(),match)!=match_identities.end()) Py_RETURN_NONE;
                    match_identities.push_back(match);
                    Owned fragment(peak_or_match_field(match,state,2,fragment_offset));
                    Owned peak(peak_or_match_field(match,state,3,peak_offset));
                    if (!fragment.value || !peak.value) return nullptr;
                    if (reinterpret_cast<PyObject*>(Py_TYPE(fragment.value))!=state->fragment_type ||
                            reinterpret_cast<PyObject*>(Py_TYPE(peak.value))!=state->peak_type) Py_RETURN_NONE;
                    const auto position=frame->position(peak.value);if (position<0) Py_RETURN_NONE;
                    const auto identity=reinterpret_cast<std::uintptr_t>(fragment.value);
                    double relative=frame->relative[position];
                    if (normalize==Py_True && frame->denominator>0.0 &&
                            std::abs(frame->masses[position]-frame->precursor)>frame->exclusion) {
                        Owned kind(model_field(fragment.value,state->fragment_fields[2],fragment_offsets[2]));
                        if (!kind.value) return nullptr;if (!PyUnicode_CheckExact(kind.value)) Py_RETURN_NONE;
                        if (PyUnicode_CompareWithASCIIString(kind.value,"Precursor Ion")!=0)
                            relative=std::max(relative,std::min(frame->raw[position]/frame->denominator,1.0));
                    }
                    candidate->matches.push_back({0,-1,reinterpret_cast<std::uintptr_t>(peak.value),identity,
                        frame->masses[position],relative,chain_flags[j],0});
                }
            } else if (prepare_arguments!=Py_None) Py_RETURN_NONE;
            remaining.push_back(i);candidates.push_back(std::move(candidate));
        }
        auto ppm_score=[&](const auto& candidate) {
            return std::max(0.0,1.0-std::min(std::abs(candidate.ppm),ppm_limit)/ppm_limit);
        };
        auto better=[&](Py_ssize_t a,Py_ssize_t b) {
            const auto& left=*candidates[a];const auto& right=*candidates[b];
            const double x=std::max(0.0,std::min(left.score/100.0,1.0));
            const double y=std::max(0.0,std::min(right.score/100.0,1.0));
            if (x!=y) return x>y;
            if (left.match_count!=right.match_count) return left.match_count>right.match_count;
            const auto p=ppm_score(left),q=ppm_score(right);if (p!=q) return p>q;
            return left.score>right.score;
        };
        while (!remaining.empty()) {
            for (auto index:remaining) {
                auto& candidate=*candidates[index];bool changed=false;
                for (const auto& group:groups) {
                    SharedPoolMatch* strongest=nullptr;
                    for (auto& match:candidate.matches) {
                        if (match.chain && std::find(group.begin(),group.end(),match.peak)!=group.end() &&
                                (!strongest || match.relative>strongest->relative)) strongest=&match;
                    }
                    if (strongest) {++strongest->penalties;changed=true;}
                }
                if (changed) {
                    if (!candidate.plan) {
                        Owned prepared(PyObject_CallObject(prepare_pool,candidate.prepare_arguments));
                        if (!prepared.value) return nullptr;
                        if (prepared.value==Py_None) Py_RETURN_NONE;
                        if (!PyTuple_CheckExact(prepared.value) || PyTuple_GET_SIZE(prepared.value)!=2) Py_RETURN_NONE;
                        auto* capsule=PyTuple_GET_ITEM(prepared.value,0);
                        auto* configuration=PyTuple_GET_ITEM(prepared.value,1);
                        candidate.plan=static_cast<PreparedPool*>(PyCapsule_GetPointer(capsule,POOL_CAPSULE));
                        if (!candidate.plan) return nullptr;
                        if (!PyTuple_CheckExact(configuration) || PyTuple_GET_SIZE(configuration)!=candidate.plan->count) Py_RETURN_NONE;
                        candidate.capsule=Py_NewRef(capsule);candidate.configuration=Py_NewRef(configuration);
                        for (auto& match:candidate.matches) {
                            auto route=std::lower_bound(candidate.plan->routes.begin(),candidate.plan->routes.end(),match.fragment,
                                [](const auto& value,auto key){return value.identity<key;});
                            if (route==candidate.plan->routes.end() || route->identity!=match.fragment) Py_RETURN_NONE;
                            match.mask=route->mask;match.group=route->group;
                        }
                    }
                    if (!shared_score(candidate,factor,score_min,score_max,sum_function,round_function)) return nullptr;
                }
            }
            remaining.erase(std::remove_if(remaining.begin(),remaining.end(),
                [&](auto index){return candidates[index]->score<minimum;}),remaining.end());
            if (remaining.empty()) break;
            std::stable_sort(remaining.begin(),remaining.end(),better);
            const double score=candidates[remaining[0]]->score;
            std::vector<Py_ssize_t> tier,after;
            for (auto index:remaining)
                (std::abs(candidates[index]->score-score)<=1e-9 ? tier : after).push_back(index);
            remaining.swap(after);selected.insert(selected.end(),tier.begin(),tier.end());
            groups.clear();
            for (auto index:tier) {
                std::vector<std::uintptr_t> group;
                for (const auto& match:candidates[index]->matches) if (match.chain &&
                        std::find(group.begin(),group.end(),match.peak)==group.end()) group.push_back(match.peak);
                if (!group.empty()) groups.push_back(std::move(group));
            }
        }
        Owned order(PyTuple_New(selected.size())),changes(PyList_New(0));
        if (!order.value || !changes.value) return nullptr;
        for (std::size_t i=0;i<selected.size();++i) {
            auto* number=PyLong_FromSsize_t(selected[i]);if (!number) return nullptr;
            PyTuple_SET_ITEM(order.value,i,number);
        }
        for (Py_ssize_t i=0;i<count;++i) if (candidates[i]->metrics) {
            Owned row(Py_BuildValue("(nOOO)",i,candidates[i]->configuration,candidates[i]->metrics,candidates[i]->score_object));
            if (!row.value || PyList_Append(changes.value,row.value)<0) return nullptr;
        }
        return PyTuple_Pack(2,order.value,changes.value);
    } catch (const std::bad_alloc&) {return PyErr_NoMemory();}
      catch (...) {PyErr_SetString(PyExc_RuntimeError,"Native shared-chain rerank failed");return nullptr;}
}
bool plain_scalar(PyObject* value, const ModuleState* state) {
    auto* type=reinterpret_cast<PyObject*>(Py_TYPE(value));
    return PyFloat_CheckExact(value) || PyLong_CheckExact(value) || PyUnicode_CheckExact(value) ||
        PyBool_Check(value) || value==Py_None || type==state->float32_type || type==state->float64_type;
}
// 1 equal, 0 edited, -1 Python error, -2 unusual public values: use Python.
int same_scalar(PyObject* current, PyObject* original, const ModuleState* state) {
    if (current==original) return 1;
    if (!plain_scalar(current,state) || !plain_scalar(original,state)) return -2;
    return PyObject_RichCompareBool(current,original,Py_EQ);
}
PyObject* unchanged_answer(int result) {
    if (result==-1) return nullptr;
    if (result==-2) Py_RETURN_NONE;
    return PyBool_FromLong(result);
}
int dictionary_unchanged(PyObject* current,PyObject* saved,const ModuleState* state) {
    if (!PyDict_CheckExact(current) || !PyTuple_CheckExact(saved)) return -2;
    if (PyDict_Size(current)!=PyTuple_GET_SIZE(saved)) return 0;
    Py_ssize_t position=0,index=0;PyObject *key,*value;
    while (PyDict_Next(current,&position,&key,&value)) {
        auto* pair=PyTuple_GET_ITEM(saved,index++);
        if (!PyTuple_CheckExact(pair) || PyTuple_GET_SIZE(pair)!=2) return -2;
        const int a=same_scalar(key,PyTuple_GET_ITEM(pair,0),state);
        if (a!=1) return a;
        const int b=same_scalar(value,PyTuple_GET_ITEM(pair,1),state);
        if (b!=1) return b;
    }
    return 1;
}
PyObject* rule_unchanged(PyObject* module,PyObject* args) {
    PyObject *rule,*revision;
    if (!PyArg_ParseTuple(args,"OO",&rule,&revision)) return nullptr;
    auto* state=static_cast<ModuleState*>(PyModule_GetState(module));
    if (reinterpret_cast<PyObject*>(Py_TYPE(rule))!=state->rule_type ||
            !PyTuple_CheckExact(revision) || PyTuple_GET_SIZE(revision)!=19) Py_RETURN_NONE;
    for (int i=0;i<15;++i) {
        Owned value(PyObject_GetAttr(rule,state->rule_fields[i]));if (!value.value) return nullptr;
        const int equal=i==1 ? dictionary_unchanged(value.value,PyTuple_GET_ITEM(revision,i),state)
                            : same_scalar(value.value,PyTuple_GET_ITEM(revision,i),state);
        if (equal!=1) return unchanged_answer(equal);
    }
    Owned profile(PyObject_GetAttr(rule,state->rule_fields[15]));if (!profile.value) return nullptr;
    if (reinterpret_cast<PyObject*>(Py_TYPE(profile.value))!=state->profile_type) Py_RETURN_NONE;
    for (int i=0;i<4;++i) {
        Owned value(PyObject_GetAttr(profile.value,state->profile_fields[i]));if (!value.value) return nullptr;
        const int equal=i<2 ? dictionary_unchanged(value.value,PyTuple_GET_ITEM(revision,15+i),state)
                            : same_scalar(value.value,PyTuple_GET_ITEM(revision,15+i),state);
        if (equal!=1) return unchanged_answer(equal);
    }
    Py_RETURN_TRUE;
}
int nested_rule_unchanged(PyObject* current,PyObject* saved,const ModuleState* state) {
    if (current==saved) return 1;
    if (PySet_CheckExact(current) || PyFrozenSet_CheckExact(current)) {
        if (!PyFrozenSet_CheckExact(saved)) return -2;
        for (auto* object:{current,saved}) {
            Owned iterator(PyObject_GetIter(object));if (!iterator.value) return -1;
            while (auto* item=PyIter_Next(iterator.value)) {
                Owned value(item);if (!PyUnicode_CheckExact(item)) return -2;
            }
            if (PyErr_Occurred()) return -1;
        }
        return PyObject_RichCompareBool(current,saved,Py_EQ);
    }
    if (PyList_CheckExact(current) || PyTuple_CheckExact(current)) {
        if (!PyTuple_CheckExact(saved)) return -2;
        if (PySequence_Fast_GET_SIZE(current)!=PyTuple_GET_SIZE(saved)) return 0;
        for (Py_ssize_t i=0;i<PyTuple_GET_SIZE(saved);++i) {
            const int equal=nested_rule_unchanged(PySequence_Fast_GET_ITEM(current,i),PyTuple_GET_ITEM(saved,i),state);
            if (equal!=1) return equal;
        }
        return 1;
    }
    if (PyDict_CheckExact(current)) {
        if (!PyTuple_CheckExact(saved)) return -2;
        if (PyDict_Size(current)!=PyTuple_GET_SIZE(saved)) return 0;
        Py_ssize_t position=0,index=0;PyObject *key,*value;
        while (PyDict_Next(current,&position,&key,&value)) {
            auto* pair=PyTuple_GET_ITEM(saved,index++);
            if (!PyTuple_CheckExact(pair) || PyTuple_GET_SIZE(pair)!=2) return -2;
            const int a=same_scalar(key,PyTuple_GET_ITEM(pair,0),state);
            if (a!=1) return a;
            const int b=nested_rule_unchanged(value,PyTuple_GET_ITEM(pair,1),state);
            if (b!=1) return b;
        }
        return 1;
    }
    return same_scalar(current,saved,state);
}
PyObject* sphingo_rule_unchanged(PyObject* module,PyObject* args) {
    PyObject *rule,*revision;
    if (!PyArg_ParseTuple(args,"OO",&rule,&revision)) return nullptr;
    auto* state=static_cast<ModuleState*>(PyModule_GetState(module));
    if (reinterpret_cast<PyObject*>(Py_TYPE(rule))!=state->sphingo_rule_type ||
            !PyTuple_CheckExact(revision) || PyTuple_GET_SIZE(revision)!=14) Py_RETURN_NONE;
    for (int i=0;i<14;++i) {
        Owned value(PyObject_GetAttr(rule,state->sphingo_rule_fields[i]));if (!value.value) return nullptr;
        const int equal=nested_rule_unchanged(value.value,PyTuple_GET_ITEM(revision,i),state);
        if (equal!=1) return unchanged_answer(equal);
    }
    Py_RETURN_TRUE;
}
PyObject* record_unchanged(PyObject* module, PyObject* args) {
    PyObject *record,*revision;
    if (!PyArg_ParseTuple(args,"OO",&record,&revision)) return nullptr;
    auto* state=static_cast<ModuleState*>(PyModule_GetState(module));
    if (reinterpret_cast<PyObject*>(Py_TYPE(record))!=state->record_type ||
            !PyTuple_CheckExact(revision) || PyTuple_GET_SIZE(revision)!=11) Py_RETURN_NONE;
    Py_ssize_t record_offsets[10], fragment_offsets[6];
    active_offsets(state,state->record_type,state->record_fields,state->record_descriptors,state->record_offsets,record_offsets);
    active_offsets(state,state->fragment_type,state->fragment_fields,state->fragment_descriptors,state->fragment_offsets,fragment_offsets);
    constexpr int field_indexes[]={0,1,2,3,5,6,7,8};
    for (int i=0; i<8; ++i) {
        Owned value(model_field(record,state->record_fields[i],record_offsets[i]));
        if (!value.value) return nullptr;
        if (i==4 && reinterpret_cast<PyObject*>(Py_TYPE(value.value))!=PyTuple_GET_ITEM(revision,4)) Py_RETURN_FALSE;
        const int equal=same_scalar(value.value,PyTuple_GET_ITEM(revision,field_indexes[i]),state);
        if (equal!=1) return unchanged_answer(equal);
    }
    Owned fragments(model_field(record,state->record_fields[8],record_offsets[8]));
    if (!fragments.value) return nullptr;
    auto* prior=PyTuple_GET_ITEM(revision,9);
    if ((!PyList_CheckExact(fragments.value) && !PyTuple_CheckExact(fragments.value)) ||
            !PyTuple_CheckExact(prior)) Py_RETURN_NONE;
    const auto count=PySequence_Fast_GET_SIZE(fragments.value);
    if (count!=PyTuple_GET_SIZE(prior)) Py_RETURN_FALSE;
    auto** current=PySequence_Fast_ITEMS(fragments.value);
    constexpr int fragment_indexes[]={2,3,4,6,8,9};
    constexpr int fragment_type_indexes[]={1,-1,-1,5,7,-1};
    for (Py_ssize_t i=0; i<count; ++i) {
        auto* saved=PyTuple_GET_ITEM(prior,i);
        if (reinterpret_cast<PyObject*>(Py_TYPE(current[i]))!=state->fragment_type ||
                !PyTuple_CheckExact(saved) || PyTuple_GET_SIZE(saved)!=10 ||
                !PyLong_CheckExact(PyTuple_GET_ITEM(saved,0))) Py_RETURN_NONE;
        auto* identity=PyLong_AsVoidPtr(PyTuple_GET_ITEM(saved,0));
        if (PyErr_Occurred()) return nullptr;
        if (identity!=current[i]) Py_RETURN_FALSE;
        for (int j=0; j<6; ++j) {
            Owned value(model_field(current[i],state->fragment_fields[j],fragment_offsets[j]));
            if (!value.value) return nullptr;
            const auto type_index=fragment_type_indexes[j];
            if (type_index>=0 && reinterpret_cast<PyObject*>(Py_TYPE(value.value))!=PyTuple_GET_ITEM(saved,type_index)) Py_RETURN_FALSE;
            const int equal=same_scalar(value.value,PyTuple_GET_ITEM(saved,fragment_indexes[j]),state);
            if (equal!=1) return unchanged_answer(equal);
        }
    }
    Owned metadata(model_field(record,state->record_fields[9],record_offsets[9]));
    if (!metadata.value) return nullptr;
    auto* saved_metadata=PyTuple_GET_ITEM(revision,10);
    if (!PyDict_CheckExact(metadata.value) || !PyTuple_CheckExact(saved_metadata)) Py_RETURN_NONE;
    if (PyDict_Size(metadata.value)!=PyTuple_GET_SIZE(saved_metadata)) Py_RETURN_FALSE;
    Py_ssize_t position=0,index=0;PyObject *key,*value;
    while (PyDict_Next(metadata.value,&position,&key,&value)) {
        auto* pair=PyTuple_GET_ITEM(saved_metadata,index++);
        if (!PyTuple_CheckExact(pair) || PyTuple_GET_SIZE(pair)!=2) Py_RETURN_NONE;
        int equal=same_scalar(key,PyTuple_GET_ITEM(pair,0),state);
        if (equal!=1) return unchanged_answer(equal);
        equal=same_scalar(value,PyTuple_GET_ITEM(pair,1),state);
        if (equal!=1) return unchanged_answer(equal);
    }
    Py_RETURN_TRUE;
}
int traverse(PyObject* module, visitproc visit, void* arg) {
    auto* s=static_cast<ModuleState*>(PyModule_GetState(module));
    Py_VISIT(s->peak_type);Py_VISIT(s->float32_type);Py_VISIT(s->float64_type);
    Py_VISIT(s->match_type);Py_VISIT(s->fragment_type);Py_VISIT(s->record_type);
    Py_VISIT(s->rule_type);Py_VISIT(s->profile_type);Py_VISIT(s->spectrum_type);Py_VISIT(s->match_initializer);
    Py_VISIT(s->sphingo_rule_type);
    for (auto* name:s->sphingo_rule_fields) Py_VISIT(name);
    for (auto* name:s->rule_fields) Py_VISIT(name);
    for (auto* name:s->profile_fields) Py_VISIT(name);
    for (auto* name:s->record_fields) Py_VISIT(name);
    for (auto* name:s->fragment_fields) Py_VISIT(name);
    for (auto* descriptor:s->record_descriptors) Py_VISIT(descriptor);
    for (auto* descriptor:s->fragment_descriptors) Py_VISIT(descriptor);
    for (auto* name:s->fast_names) Py_VISIT(name);
    for (auto* descriptor:s->fast_descriptors) Py_VISIT(descriptor);
    return 0;
}
int clear(PyObject* module) {
    auto* s=static_cast<ModuleState*>(PyModule_GetState(module));
    Py_CLEAR(s->peak_type);Py_CLEAR(s->float32_type);Py_CLEAR(s->float64_type);
    Py_CLEAR(s->match_type);Py_CLEAR(s->fragment_type);Py_CLEAR(s->record_type);
    Py_CLEAR(s->rule_type);Py_CLEAR(s->profile_type);Py_CLEAR(s->spectrum_type);Py_CLEAR(s->match_initializer);
    Py_CLEAR(s->sphingo_rule_type);
    for (auto*& name:s->sphingo_rule_fields) Py_CLEAR(name);
    for (auto*& name:s->rule_fields) Py_CLEAR(name);
    for (auto*& name:s->profile_fields) Py_CLEAR(name);
    for (auto*& name:s->record_fields) Py_CLEAR(name);
    for (auto*& name:s->fragment_fields) Py_CLEAR(name);
    for (auto*& descriptor:s->record_descriptors) Py_CLEAR(descriptor);
    for (auto*& descriptor:s->fragment_descriptors) Py_CLEAR(descriptor);
    for (auto*& name:s->fast_names) Py_CLEAR(name);
    for (auto*& descriptor:s->fast_descriptors) Py_CLEAR(descriptor);
    return 0;
}
PyMethodDef methods[] = {{"create_plan",create_plan,METH_VARARGS,nullptr},
                        {"match_plan",match_plan,METH_VARARGS,nullptr},
                        {"match_query_plan",match_query_plan,METH_VARARGS,nullptr},
                        {"pool_metrics",pool_metrics,METH_VARARGS,nullptr},
                        {"pool_metrics_query",pool_metrics,METH_VARARGS,nullptr},
                        {"score_pool_policy",score_pool_policy,METH_VARARGS,nullptr},
                        {"pool_metrics_normalized",pool_metrics_normalized,METH_VARARGS,nullptr},
                        {"score_pool_policy_normalized",score_pool_policy_normalized,METH_VARARGS,nullptr},
                        {"score_pool_indexes",score_pool_indexes,METH_VARARGS,nullptr},
                        {"pool_metrics_indexes",pool_metrics_indexes,METH_VARARGS,nullptr},
                        {"indexed_scoring_compatible",indexed_scoring_compatible,METH_O,nullptr},
                        {"spectrum_is_builtin",spectrum_is_builtin,METH_O,nullptr},
                        {"non_precursor_base",non_precursor_base,METH_VARARGS,nullptr},
                        {"materialize_matches",materialize_matches,METH_VARARGS,nullptr},
                        {"materialize_query_matches",materialize_query_matches,METH_VARARGS,nullptr},
                        {"quality_overrides_query",quality_overrides_query,METH_VARARGS,nullptr},
                        {"query_base",query_base,METH_O,nullptr},
                        {"query_intensity_sums",query_intensity_sums,METH_VARARGS,nullptr},
                        {"create_pool_plan",create_pool_plan,METH_VARARGS,nullptr},
                        {"quality_overrides",quality_overrides,METH_VARARGS,nullptr},
                        {"rule_unchanged",rule_unchanged,METH_VARARGS,nullptr},
                        {"sphingo_rule_unchanged",sphingo_rule_unchanged,METH_VARARGS,nullptr},
                        {"record_unchanged",record_unchanged,METH_VARARGS,nullptr},
                        {"rerank_pool_candidates",rerank_pool_candidates,METH_VARARGS,nullptr},
                        {nullptr,nullptr,0,nullptr}};
PyModuleDef definition = {PyModuleDef_HEAD_INIT,"_search_bridge",nullptr,sizeof(ModuleState),methods,
                         nullptr,traverse,clear,nullptr};
}
PyMODINIT_FUNC PyInit__search_bridge() {
    PyObject* module = PyModule_Create(&definition);if (!module) return nullptr;
    auto* s=static_cast<ModuleState*>(PyModule_GetState(module));
    PyObject* models=PyImport_ImportModule("lipidgate.ms2.models");
    if (models) {
        s->peak_type=PyObject_GetAttrString(models,"ExperimentalPeak");
        s->match_type=PyObject_GetAttrString(models,"FragmentMatch");
        s->fragment_type=PyObject_GetAttrString(models,"FragmentRecord");
        s->record_type=PyObject_GetAttrString(models,"LibraryRecord");
        s->spectrum_type=PyObject_GetAttrString(models,"ExperimentalSpectrum");
        s->match_initializer=PyObject_GetAttrString(models,"_FRAGMENT_MATCH_INITIALIZER");
        if (!s->match_initializer && PyErr_ExceptionMatches(PyExc_AttributeError)) PyErr_Clear();
        Py_DECREF(models);
    }
    PyObject* numpy=PyImport_ImportModule("numpy");
    if (numpy) {
        s->float32_type=PyObject_GetAttrString(numpy,"float32");
        s->float64_type=PyObject_GetAttrString(numpy,"float64");Py_DECREF(numpy);
    }
    if (!s->peak_type || !s->float32_type || !s->float64_type || !s->match_type || !s->fragment_type || !s->record_type || !s->spectrum_type) { Py_DECREF(module);return nullptr; }
    Owned rules(PyImport_ImportModule("lipidgate.ms2.rules"));if (!rules.value) {Py_DECREF(module);return nullptr;}
    s->rule_type=PyObject_GetAttrString(rules.value,"ClassRule");
    s->profile_type=PyObject_GetAttrString(rules.value,"ScoreProfile");
    if (!s->rule_type || !s->profile_type) {Py_DECREF(module);return nullptr;}
    Owned sphingo_rules(PyImport_ImportModule("lipidgate.ms2.sphingolipid_rules"));
    if (!sphingo_rules.value) {Py_DECREF(module);return nullptr;}
    s->sphingo_rule_type=PyObject_GetAttrString(sphingo_rules.value,"SphingoRule");
    if (!s->sphingo_rule_type) {Py_DECREF(module);return nullptr;}
    for (int i=0;i<14;++i) {
        s->sphingo_rule_fields[i]=PyUnicode_InternFromString(SPHINGO_RULE_FIELDS[i]);
        if (!s->sphingo_rule_fields[i]) {Py_DECREF(module);return nullptr;}
    }
    for (int i=0;i<16;++i) {
        s->rule_fields[i]=PyUnicode_InternFromString(RULE_FIELDS[i]);
        if (!s->rule_fields[i]) {Py_DECREF(module);return nullptr;}
    }
    for (int i=0;i<4;++i) {
        s->profile_fields[i]=PyUnicode_InternFromString(PROFILE_FIELDS[i]);
        if (!s->profile_fields[i]) {Py_DECREF(module);return nullptr;}
    }
    for (int i=0; i<10; ++i) {
        s->record_fields[i]=PyUnicode_InternFromString(RECORD_FIELDS[i]);
        if (!s->record_fields[i]) {Py_DECREF(module);return nullptr;}
        auto* type=reinterpret_cast<PyTypeObject*>(s->record_type);
        auto* descriptor=PyDict_GetItemWithError(type->tp_dict,s->record_fields[i]);
        if (PyErr_Occurred()) {Py_DECREF(module);return nullptr;}
        s->record_descriptors[i]=Py_XNewRef(descriptor);
        s->record_offsets[i]=descriptor_offset(type,descriptor);
    }
    for (int i=0; i<6; ++i) {
        s->fragment_fields[i]=PyUnicode_InternFromString(FRAGMENT_FIELDS[i]);
        if (!s->fragment_fields[i]) {Py_DECREF(module);return nullptr;}
        auto* type=reinterpret_cast<PyTypeObject*>(s->fragment_type);
        auto* descriptor=PyDict_GetItemWithError(type->tp_dict,s->fragment_fields[i]);
        if (PyErr_Occurred()) {Py_DECREF(module);return nullptr;}
        s->fragment_descriptors[i]=Py_XNewRef(descriptor);
        s->fragment_offsets[i]=descriptor_offset(type,descriptor);
    }
    for (int i=0; i<6; ++i) {
        s->fast_offsets[i]=-1;
        s->fast_names[i]=PyUnicode_InternFromString(FAST_FIELDS[i]);
        if (!s->fast_names[i]) {Py_DECREF(module);return nullptr;}
        auto* type=reinterpret_cast<PyTypeObject*>(i<2 || i==4 ? s->peak_type : s->match_type);
        auto* descriptor=PyDict_GetItemWithError(type->tp_dict,s->fast_names[i]);
        if (descriptor) {
            s->fast_descriptors[i]=Py_NewRef(descriptor);
            s->fast_offsets[i]=descriptor_offset(type,descriptor);
        } else if (PyErr_Occurred()) {Py_DECREF(module);return nullptr;}
    }
    return module;
}
