// Generated by the protocol buffer compiler.  DO NOT EDIT!

#define INTERNAL_SUPPRESS_PROTOBUF_FIELD_DEPRECATION
#include "indexpb.h"

#include <algorithm>

#include <google/protobuf/stubs/once.h>
#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/wire_format_lite_inl.h>
#include <google/protobuf/descriptor.h>
#include <google/protobuf/reflection_ops.h>
#include <google/protobuf/wire_format.h>
// @@protoc_insertion_point(includes)

namespace idxfile {

namespace {

const ::google::protobuf::Descriptor* SigUnit_descriptor_ = NULL;
const ::google::protobuf::internal::GeneratedMessageReflection*
  SigUnit_reflection_ = NULL;
const ::google::protobuf::Descriptor* Entry_descriptor_ = NULL;
const ::google::protobuf::internal::GeneratedMessageReflection*
  Entry_reflection_ = NULL;
const ::google::protobuf::Descriptor* EntryList_descriptor_ = NULL;
const ::google::protobuf::internal::GeneratedMessageReflection*
  EntryList_reflection_ = NULL;

}  // namespace


void protobuf_AssignDesc_index_2eproto() {
  protobuf_AddDesc_index_2eproto();
  const ::google::protobuf::FileDescriptor* file =
    ::google::protobuf::DescriptorPool::generated_pool()->FindFileByName(
      "index.proto");
  GOOGLE_CHECK(file != NULL);
  SigUnit_descriptor_ = file->message_type(0);
  static const int SigUnit_offsets_[3] = {
    GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(SigUnit, init_),
    GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(SigUnit, deltas_),
    GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(SigUnit, cnt_),
  };
  SigUnit_reflection_ =
    new ::google::protobuf::internal::GeneratedMessageReflection(
      SigUnit_descriptor_,
      SigUnit::default_instance_,
      SigUnit_offsets_,
      GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(SigUnit, _has_bits_[0]),
      GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(SigUnit, _unknown_fields_),
      -1,
      ::google::protobuf::DescriptorPool::generated_pool(),
      ::google::protobuf::MessageFactory::generated_factory(),
      sizeof(SigUnit));
  Entry_descriptor_ = file->message_type(1);
  static const int Entry_offsets_[4] = {
    GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(Entry, proc_),
    GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(Entry, logical_offset_),
    GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(Entry, length_),
    GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(Entry, physical_offset_),
  };
  Entry_reflection_ =
    new ::google::protobuf::internal::GeneratedMessageReflection(
      Entry_descriptor_,
      Entry::default_instance_,
      Entry_offsets_,
      GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(Entry, _has_bits_[0]),
      GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(Entry, _unknown_fields_),
      -1,
      ::google::protobuf::DescriptorPool::generated_pool(),
      ::google::protobuf::MessageFactory::generated_factory(),
      sizeof(Entry));
  EntryList_descriptor_ = file->message_type(2);
  static const int EntryList_offsets_[1] = {
    GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(EntryList, entry_),
  };
  EntryList_reflection_ =
    new ::google::protobuf::internal::GeneratedMessageReflection(
      EntryList_descriptor_,
      EntryList::default_instance_,
      EntryList_offsets_,
      GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(EntryList, _has_bits_[0]),
      GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(EntryList, _unknown_fields_),
      -1,
      ::google::protobuf::DescriptorPool::generated_pool(),
      ::google::protobuf::MessageFactory::generated_factory(),
      sizeof(EntryList));
}

namespace {

GOOGLE_PROTOBUF_DECLARE_ONCE(protobuf_AssignDescriptors_once_);
inline void protobuf_AssignDescriptorsOnce() {
  ::google::protobuf::GoogleOnceInit(&protobuf_AssignDescriptors_once_,
                 &protobuf_AssignDesc_index_2eproto);
}

void protobuf_RegisterTypes(const ::std::string&) {
  protobuf_AssignDescriptorsOnce();
  ::google::protobuf::MessageFactory::InternalRegisterGeneratedMessage(
    SigUnit_descriptor_, &SigUnit::default_instance());
  ::google::protobuf::MessageFactory::InternalRegisterGeneratedMessage(
    Entry_descriptor_, &Entry::default_instance());
  ::google::protobuf::MessageFactory::InternalRegisterGeneratedMessage(
    EntryList_descriptor_, &EntryList::default_instance());
}

}  // namespace

void protobuf_ShutdownFile_index_2eproto() {
  delete SigUnit::default_instance_;
  delete SigUnit_reflection_;
  delete Entry::default_instance_;
  delete Entry_reflection_;
  delete EntryList::default_instance_;
  delete EntryList_reflection_;
}

void protobuf_AddDesc_index_2eproto() {
  static bool already_here = false;
  if (already_here) return;
  already_here = true;
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  ::google::protobuf::DescriptorPool::InternalAddGeneratedFile(
    "\n\013index.proto\022\007idxfile\"4\n\007SigUnit\022\014\n\004ini"
    "t\030\001 \002(\003\022\016\n\006deltas\030\002 \003(\003\022\013\n\003cnt\030\003 \002(\003\"\214\001\n"
    "\005Entry\022\014\n\004proc\030\001 \002(\003\022(\n\016logical_offset\030\002"
    " \002(\0132\020.idxfile.SigUnit\022 \n\006length\030\003 \003(\0132\020"
    ".idxfile.SigUnit\022)\n\017physical_offset\030\004 \003("
    "\0132\020.idxfile.SigUnit\"*\n\tEntryList\022\035\n\005entr"
    "y\030\001 \003(\0132\016.idxfile.Entry", 263);
  ::google::protobuf::MessageFactory::InternalRegisterGeneratedFile(
    "index.proto", &protobuf_RegisterTypes);
  SigUnit::default_instance_ = new SigUnit();
  Entry::default_instance_ = new Entry();
  EntryList::default_instance_ = new EntryList();
  SigUnit::default_instance_->InitAsDefaultInstance();
  Entry::default_instance_->InitAsDefaultInstance();
  EntryList::default_instance_->InitAsDefaultInstance();
  ::google::protobuf::internal::OnShutdown(&protobuf_ShutdownFile_index_2eproto);
}

// Force AddDescriptors() to be called at static initialization time.
struct StaticDescriptorInitializer_index_2eproto {
  StaticDescriptorInitializer_index_2eproto() {
    protobuf_AddDesc_index_2eproto();
  }
} static_descriptor_initializer_index_2eproto_;


// ===================================================================

#ifndef _MSC_VER
const int SigUnit::kInitFieldNumber;
const int SigUnit::kDeltasFieldNumber;
const int SigUnit::kCntFieldNumber;
#endif  // !_MSC_VER

SigUnit::SigUnit()
  : ::google::protobuf::Message() {
  SharedCtor();
}

void SigUnit::InitAsDefaultInstance() {
}

SigUnit::SigUnit(const SigUnit& from)
  : ::google::protobuf::Message() {
  SharedCtor();
  MergeFrom(from);
}

void SigUnit::SharedCtor() {
  _cached_size_ = 0;
  init_ = GOOGLE_LONGLONG(0);
  cnt_ = GOOGLE_LONGLONG(0);
  ::memset(_has_bits_, 0, sizeof(_has_bits_));
}

SigUnit::~SigUnit() {
  SharedDtor();
}

void SigUnit::SharedDtor() {
  if (this != default_instance_) {
  }
}

void SigUnit::SetCachedSize(int size) const {
  GOOGLE_SAFE_CONCURRENT_WRITES_BEGIN();
  _cached_size_ = size;
  GOOGLE_SAFE_CONCURRENT_WRITES_END();
}
const ::google::protobuf::Descriptor* SigUnit::descriptor() {
  protobuf_AssignDescriptorsOnce();
  return SigUnit_descriptor_;
}

const SigUnit& SigUnit::default_instance() {
  if (default_instance_ == NULL) protobuf_AddDesc_index_2eproto();  return *default_instance_;
}

SigUnit* SigUnit::default_instance_ = NULL;

SigUnit* SigUnit::New() const {
  return new SigUnit;
}

void SigUnit::Clear() {
  if (_has_bits_[0 / 32] & (0xffu << (0 % 32))) {
    init_ = GOOGLE_LONGLONG(0);
    cnt_ = GOOGLE_LONGLONG(0);
  }
  deltas_.Clear();
  ::memset(_has_bits_, 0, sizeof(_has_bits_));
  mutable_unknown_fields()->Clear();
}

bool SigUnit::MergePartialFromCodedStream(
    ::google::protobuf::io::CodedInputStream* input) {
#define DO_(EXPRESSION) if (!(EXPRESSION)) return false
  ::google::protobuf::uint32 tag;
  while ((tag = input->ReadTag()) != 0) {
    switch (::google::protobuf::internal::WireFormatLite::GetTagFieldNumber(tag)) {
      // required int64 init = 1;
      case 1: {
        if (::google::protobuf::internal::WireFormatLite::GetTagWireType(tag) ==
            ::google::protobuf::internal::WireFormatLite::WIRETYPE_VARINT) {
          DO_((::google::protobuf::internal::WireFormatLite::ReadPrimitive<
                   ::google::protobuf::int64, ::google::protobuf::internal::WireFormatLite::TYPE_INT64>(
                 input, &init_)));
          set_has_init();
        } else {
          goto handle_uninterpreted;
        }
        if (input->ExpectTag(16)) goto parse_deltas;
        break;
      }
      
      // repeated int64 deltas = 2;
      case 2: {
        if (::google::protobuf::internal::WireFormatLite::GetTagWireType(tag) ==
            ::google::protobuf::internal::WireFormatLite::WIRETYPE_VARINT) {
         parse_deltas:
          DO_((::google::protobuf::internal::WireFormatLite::ReadRepeatedPrimitive<
                   ::google::protobuf::int64, ::google::protobuf::internal::WireFormatLite::TYPE_INT64>(
                 1, 16, input, this->mutable_deltas())));
        } else if (::google::protobuf::internal::WireFormatLite::GetTagWireType(tag)
                   == ::google::protobuf::internal::WireFormatLite::
                      WIRETYPE_LENGTH_DELIMITED) {
          DO_((::google::protobuf::internal::WireFormatLite::ReadPackedPrimitiveNoInline<
                   ::google::protobuf::int64, ::google::protobuf::internal::WireFormatLite::TYPE_INT64>(
                 input, this->mutable_deltas())));
        } else {
          goto handle_uninterpreted;
        }
        if (input->ExpectTag(16)) goto parse_deltas;
        if (input->ExpectTag(24)) goto parse_cnt;
        break;
      }
      
      // required int64 cnt = 3;
      case 3: {
        if (::google::protobuf::internal::WireFormatLite::GetTagWireType(tag) ==
            ::google::protobuf::internal::WireFormatLite::WIRETYPE_VARINT) {
         parse_cnt:
          DO_((::google::protobuf::internal::WireFormatLite::ReadPrimitive<
                   ::google::protobuf::int64, ::google::protobuf::internal::WireFormatLite::TYPE_INT64>(
                 input, &cnt_)));
          set_has_cnt();
        } else {
          goto handle_uninterpreted;
        }
        if (input->ExpectAtEnd()) return true;
        break;
      }
      
      default: {
      handle_uninterpreted:
        if (::google::protobuf::internal::WireFormatLite::GetTagWireType(tag) ==
            ::google::protobuf::internal::WireFormatLite::WIRETYPE_END_GROUP) {
          return true;
        }
        DO_(::google::protobuf::internal::WireFormat::SkipField(
              input, tag, mutable_unknown_fields()));
        break;
      }
    }
  }
  return true;
#undef DO_
}

void SigUnit::SerializeWithCachedSizes(
    ::google::protobuf::io::CodedOutputStream* output) const {
  // required int64 init = 1;
  if (has_init()) {
    ::google::protobuf::internal::WireFormatLite::WriteInt64(1, this->init(), output);
  }
  
  // repeated int64 deltas = 2;
  for (int i = 0; i < this->deltas_size(); i++) {
    ::google::protobuf::internal::WireFormatLite::WriteInt64(
      2, this->deltas(i), output);
  }
  
  // required int64 cnt = 3;
  if (has_cnt()) {
    ::google::protobuf::internal::WireFormatLite::WriteInt64(3, this->cnt(), output);
  }
  
  if (!unknown_fields().empty()) {
    ::google::protobuf::internal::WireFormat::SerializeUnknownFields(
        unknown_fields(), output);
  }
}

::google::protobuf::uint8* SigUnit::SerializeWithCachedSizesToArray(
    ::google::protobuf::uint8* target) const {
  // required int64 init = 1;
  if (has_init()) {
    target = ::google::protobuf::internal::WireFormatLite::WriteInt64ToArray(1, this->init(), target);
  }
  
  // repeated int64 deltas = 2;
  for (int i = 0; i < this->deltas_size(); i++) {
    target = ::google::protobuf::internal::WireFormatLite::
      WriteInt64ToArray(2, this->deltas(i), target);
  }
  
  // required int64 cnt = 3;
  if (has_cnt()) {
    target = ::google::protobuf::internal::WireFormatLite::WriteInt64ToArray(3, this->cnt(), target);
  }
  
  if (!unknown_fields().empty()) {
    target = ::google::protobuf::internal::WireFormat::SerializeUnknownFieldsToArray(
        unknown_fields(), target);
  }
  return target;
}

int SigUnit::ByteSize() const {
  int total_size = 0;
  
  if (_has_bits_[0 / 32] & (0xffu << (0 % 32))) {
    // required int64 init = 1;
    if (has_init()) {
      total_size += 1 +
        ::google::protobuf::internal::WireFormatLite::Int64Size(
          this->init());
    }
    
    // required int64 cnt = 3;
    if (has_cnt()) {
      total_size += 1 +
        ::google::protobuf::internal::WireFormatLite::Int64Size(
          this->cnt());
    }
    
  }
  // repeated int64 deltas = 2;
  {
    int data_size = 0;
    for (int i = 0; i < this->deltas_size(); i++) {
      data_size += ::google::protobuf::internal::WireFormatLite::
        Int64Size(this->deltas(i));
    }
    total_size += 1 * this->deltas_size() + data_size;
  }
  
  if (!unknown_fields().empty()) {
    total_size +=
      ::google::protobuf::internal::WireFormat::ComputeUnknownFieldsSize(
        unknown_fields());
  }
  GOOGLE_SAFE_CONCURRENT_WRITES_BEGIN();
  _cached_size_ = total_size;
  GOOGLE_SAFE_CONCURRENT_WRITES_END();
  return total_size;
}

void SigUnit::MergeFrom(const ::google::protobuf::Message& from) {
  GOOGLE_CHECK_NE(&from, this);
  const SigUnit* source =
    ::google::protobuf::internal::dynamic_cast_if_available<const SigUnit*>(
      &from);
  if (source == NULL) {
    ::google::protobuf::internal::ReflectionOps::Merge(from, this);
  } else {
    MergeFrom(*source);
  }
}

void SigUnit::MergeFrom(const SigUnit& from) {
  GOOGLE_CHECK_NE(&from, this);
  deltas_.MergeFrom(from.deltas_);
  if (from._has_bits_[0 / 32] & (0xffu << (0 % 32))) {
    if (from.has_init()) {
      set_init(from.init());
    }
    if (from.has_cnt()) {
      set_cnt(from.cnt());
    }
  }
  mutable_unknown_fields()->MergeFrom(from.unknown_fields());
}

void SigUnit::CopyFrom(const ::google::protobuf::Message& from) {
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

void SigUnit::CopyFrom(const SigUnit& from) {
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool SigUnit::IsInitialized() const {
  if ((_has_bits_[0] & 0x00000005) != 0x00000005) return false;
  
  return true;
}

void SigUnit::Swap(SigUnit* other) {
  if (other != this) {
    std::swap(init_, other->init_);
    deltas_.Swap(&other->deltas_);
    std::swap(cnt_, other->cnt_);
    std::swap(_has_bits_[0], other->_has_bits_[0]);
    _unknown_fields_.Swap(&other->_unknown_fields_);
    std::swap(_cached_size_, other->_cached_size_);
  }
}

::google::protobuf::Metadata SigUnit::GetMetadata() const {
  protobuf_AssignDescriptorsOnce();
  ::google::protobuf::Metadata metadata;
  metadata.descriptor = SigUnit_descriptor_;
  metadata.reflection = SigUnit_reflection_;
  return metadata;
}


// ===================================================================

#ifndef _MSC_VER
const int Entry::kProcFieldNumber;
const int Entry::kLogicalOffsetFieldNumber;
const int Entry::kLengthFieldNumber;
const int Entry::kPhysicalOffsetFieldNumber;
#endif  // !_MSC_VER

Entry::Entry()
  : ::google::protobuf::Message() {
  SharedCtor();
}

void Entry::InitAsDefaultInstance() {
  logical_offset_ = const_cast< ::idxfile::SigUnit*>(&::idxfile::SigUnit::default_instance());
}

Entry::Entry(const Entry& from)
  : ::google::protobuf::Message() {
  SharedCtor();
  MergeFrom(from);
}

void Entry::SharedCtor() {
  _cached_size_ = 0;
  proc_ = GOOGLE_LONGLONG(0);
  logical_offset_ = NULL;
  ::memset(_has_bits_, 0, sizeof(_has_bits_));
}

Entry::~Entry() {
  SharedDtor();
}

void Entry::SharedDtor() {
  if (this != default_instance_) {
    delete logical_offset_;
  }
}

void Entry::SetCachedSize(int size) const {
  GOOGLE_SAFE_CONCURRENT_WRITES_BEGIN();
  _cached_size_ = size;
  GOOGLE_SAFE_CONCURRENT_WRITES_END();
}
const ::google::protobuf::Descriptor* Entry::descriptor() {
  protobuf_AssignDescriptorsOnce();
  return Entry_descriptor_;
}

const Entry& Entry::default_instance() {
  if (default_instance_ == NULL) protobuf_AddDesc_index_2eproto();  return *default_instance_;
}

Entry* Entry::default_instance_ = NULL;

Entry* Entry::New() const {
  return new Entry;
}

void Entry::Clear() {
  if (_has_bits_[0 / 32] & (0xffu << (0 % 32))) {
    proc_ = GOOGLE_LONGLONG(0);
    if (has_logical_offset()) {
      if (logical_offset_ != NULL) logical_offset_->::idxfile::SigUnit::Clear();
    }
  }
  length_.Clear();
  physical_offset_.Clear();
  ::memset(_has_bits_, 0, sizeof(_has_bits_));
  mutable_unknown_fields()->Clear();
}

bool Entry::MergePartialFromCodedStream(
    ::google::protobuf::io::CodedInputStream* input) {
#define DO_(EXPRESSION) if (!(EXPRESSION)) return false
  ::google::protobuf::uint32 tag;
  while ((tag = input->ReadTag()) != 0) {
    switch (::google::protobuf::internal::WireFormatLite::GetTagFieldNumber(tag)) {
      // required int64 proc = 1;
      case 1: {
        if (::google::protobuf::internal::WireFormatLite::GetTagWireType(tag) ==
            ::google::protobuf::internal::WireFormatLite::WIRETYPE_VARINT) {
          DO_((::google::protobuf::internal::WireFormatLite::ReadPrimitive<
                   ::google::protobuf::int64, ::google::protobuf::internal::WireFormatLite::TYPE_INT64>(
                 input, &proc_)));
          set_has_proc();
        } else {
          goto handle_uninterpreted;
        }
        if (input->ExpectTag(18)) goto parse_logical_offset;
        break;
      }
      
      // required .idxfile.SigUnit logical_offset = 2;
      case 2: {
        if (::google::protobuf::internal::WireFormatLite::GetTagWireType(tag) ==
            ::google::protobuf::internal::WireFormatLite::WIRETYPE_LENGTH_DELIMITED) {
         parse_logical_offset:
          DO_(::google::protobuf::internal::WireFormatLite::ReadMessageNoVirtual(
               input, mutable_logical_offset()));
        } else {
          goto handle_uninterpreted;
        }
        if (input->ExpectTag(26)) goto parse_length;
        break;
      }
      
      // repeated .idxfile.SigUnit length = 3;
      case 3: {
        if (::google::protobuf::internal::WireFormatLite::GetTagWireType(tag) ==
            ::google::protobuf::internal::WireFormatLite::WIRETYPE_LENGTH_DELIMITED) {
         parse_length:
          DO_(::google::protobuf::internal::WireFormatLite::ReadMessageNoVirtual(
                input, add_length()));
        } else {
          goto handle_uninterpreted;
        }
        if (input->ExpectTag(26)) goto parse_length;
        if (input->ExpectTag(34)) goto parse_physical_offset;
        break;
      }
      
      // repeated .idxfile.SigUnit physical_offset = 4;
      case 4: {
        if (::google::protobuf::internal::WireFormatLite::GetTagWireType(tag) ==
            ::google::protobuf::internal::WireFormatLite::WIRETYPE_LENGTH_DELIMITED) {
         parse_physical_offset:
          DO_(::google::protobuf::internal::WireFormatLite::ReadMessageNoVirtual(
                input, add_physical_offset()));
        } else {
          goto handle_uninterpreted;
        }
        if (input->ExpectTag(34)) goto parse_physical_offset;
        if (input->ExpectAtEnd()) return true;
        break;
      }
      
      default: {
      handle_uninterpreted:
        if (::google::protobuf::internal::WireFormatLite::GetTagWireType(tag) ==
            ::google::protobuf::internal::WireFormatLite::WIRETYPE_END_GROUP) {
          return true;
        }
        DO_(::google::protobuf::internal::WireFormat::SkipField(
              input, tag, mutable_unknown_fields()));
        break;
      }
    }
  }
  return true;
#undef DO_
}

void Entry::SerializeWithCachedSizes(
    ::google::protobuf::io::CodedOutputStream* output) const {
  // required int64 proc = 1;
  if (has_proc()) {
    ::google::protobuf::internal::WireFormatLite::WriteInt64(1, this->proc(), output);
  }
  
  // required .idxfile.SigUnit logical_offset = 2;
  if (has_logical_offset()) {
    ::google::protobuf::internal::WireFormatLite::WriteMessageMaybeToArray(
      2, this->logical_offset(), output);
  }
  
  // repeated .idxfile.SigUnit length = 3;
  for (int i = 0; i < this->length_size(); i++) {
    ::google::protobuf::internal::WireFormatLite::WriteMessageMaybeToArray(
      3, this->length(i), output);
  }
  
  // repeated .idxfile.SigUnit physical_offset = 4;
  for (int i = 0; i < this->physical_offset_size(); i++) {
    ::google::protobuf::internal::WireFormatLite::WriteMessageMaybeToArray(
      4, this->physical_offset(i), output);
  }
  
  if (!unknown_fields().empty()) {
    ::google::protobuf::internal::WireFormat::SerializeUnknownFields(
        unknown_fields(), output);
  }
}

::google::protobuf::uint8* Entry::SerializeWithCachedSizesToArray(
    ::google::protobuf::uint8* target) const {
  // required int64 proc = 1;
  if (has_proc()) {
    target = ::google::protobuf::internal::WireFormatLite::WriteInt64ToArray(1, this->proc(), target);
  }
  
  // required .idxfile.SigUnit logical_offset = 2;
  if (has_logical_offset()) {
    target = ::google::protobuf::internal::WireFormatLite::
      WriteMessageNoVirtualToArray(
        2, this->logical_offset(), target);
  }
  
  // repeated .idxfile.SigUnit length = 3;
  for (int i = 0; i < this->length_size(); i++) {
    target = ::google::protobuf::internal::WireFormatLite::
      WriteMessageNoVirtualToArray(
        3, this->length(i), target);
  }
  
  // repeated .idxfile.SigUnit physical_offset = 4;
  for (int i = 0; i < this->physical_offset_size(); i++) {
    target = ::google::protobuf::internal::WireFormatLite::
      WriteMessageNoVirtualToArray(
        4, this->physical_offset(i), target);
  }
  
  if (!unknown_fields().empty()) {
    target = ::google::protobuf::internal::WireFormat::SerializeUnknownFieldsToArray(
        unknown_fields(), target);
  }
  return target;
}

int Entry::ByteSize() const {
  int total_size = 0;
  
  if (_has_bits_[0 / 32] & (0xffu << (0 % 32))) {
    // required int64 proc = 1;
    if (has_proc()) {
      total_size += 1 +
        ::google::protobuf::internal::WireFormatLite::Int64Size(
          this->proc());
    }
    
    // required .idxfile.SigUnit logical_offset = 2;
    if (has_logical_offset()) {
      total_size += 1 +
        ::google::protobuf::internal::WireFormatLite::MessageSizeNoVirtual(
          this->logical_offset());
    }
    
  }
  // repeated .idxfile.SigUnit length = 3;
  total_size += 1 * this->length_size();
  for (int i = 0; i < this->length_size(); i++) {
    total_size +=
      ::google::protobuf::internal::WireFormatLite::MessageSizeNoVirtual(
        this->length(i));
  }
  
  // repeated .idxfile.SigUnit physical_offset = 4;
  total_size += 1 * this->physical_offset_size();
  for (int i = 0; i < this->physical_offset_size(); i++) {
    total_size +=
      ::google::protobuf::internal::WireFormatLite::MessageSizeNoVirtual(
        this->physical_offset(i));
  }
  
  if (!unknown_fields().empty()) {
    total_size +=
      ::google::protobuf::internal::WireFormat::ComputeUnknownFieldsSize(
        unknown_fields());
  }
  GOOGLE_SAFE_CONCURRENT_WRITES_BEGIN();
  _cached_size_ = total_size;
  GOOGLE_SAFE_CONCURRENT_WRITES_END();
  return total_size;
}

void Entry::MergeFrom(const ::google::protobuf::Message& from) {
  GOOGLE_CHECK_NE(&from, this);
  const Entry* source =
    ::google::protobuf::internal::dynamic_cast_if_available<const Entry*>(
      &from);
  if (source == NULL) {
    ::google::protobuf::internal::ReflectionOps::Merge(from, this);
  } else {
    MergeFrom(*source);
  }
}

void Entry::MergeFrom(const Entry& from) {
  GOOGLE_CHECK_NE(&from, this);
  length_.MergeFrom(from.length_);
  physical_offset_.MergeFrom(from.physical_offset_);
  if (from._has_bits_[0 / 32] & (0xffu << (0 % 32))) {
    if (from.has_proc()) {
      set_proc(from.proc());
    }
    if (from.has_logical_offset()) {
      mutable_logical_offset()->::idxfile::SigUnit::MergeFrom(from.logical_offset());
    }
  }
  mutable_unknown_fields()->MergeFrom(from.unknown_fields());
}

void Entry::CopyFrom(const ::google::protobuf::Message& from) {
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

void Entry::CopyFrom(const Entry& from) {
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool Entry::IsInitialized() const {
  if ((_has_bits_[0] & 0x00000003) != 0x00000003) return false;
  
  if (has_logical_offset()) {
    if (!this->logical_offset().IsInitialized()) return false;
  }
  for (int i = 0; i < length_size(); i++) {
    if (!this->length(i).IsInitialized()) return false;
  }
  for (int i = 0; i < physical_offset_size(); i++) {
    if (!this->physical_offset(i).IsInitialized()) return false;
  }
  return true;
}

void Entry::Swap(Entry* other) {
  if (other != this) {
    std::swap(proc_, other->proc_);
    std::swap(logical_offset_, other->logical_offset_);
    length_.Swap(&other->length_);
    physical_offset_.Swap(&other->physical_offset_);
    std::swap(_has_bits_[0], other->_has_bits_[0]);
    _unknown_fields_.Swap(&other->_unknown_fields_);
    std::swap(_cached_size_, other->_cached_size_);
  }
}

::google::protobuf::Metadata Entry::GetMetadata() const {
  protobuf_AssignDescriptorsOnce();
  ::google::protobuf::Metadata metadata;
  metadata.descriptor = Entry_descriptor_;
  metadata.reflection = Entry_reflection_;
  return metadata;
}


// ===================================================================

#ifndef _MSC_VER
const int EntryList::kEntryFieldNumber;
#endif  // !_MSC_VER

EntryList::EntryList()
  : ::google::protobuf::Message() {
  SharedCtor();
}

void EntryList::InitAsDefaultInstance() {
}

EntryList::EntryList(const EntryList& from)
  : ::google::protobuf::Message() {
  SharedCtor();
  MergeFrom(from);
}

void EntryList::SharedCtor() {
  _cached_size_ = 0;
  ::memset(_has_bits_, 0, sizeof(_has_bits_));
}

EntryList::~EntryList() {
  SharedDtor();
}

void EntryList::SharedDtor() {
  if (this != default_instance_) {
  }
}

void EntryList::SetCachedSize(int size) const {
  GOOGLE_SAFE_CONCURRENT_WRITES_BEGIN();
  _cached_size_ = size;
  GOOGLE_SAFE_CONCURRENT_WRITES_END();
}
const ::google::protobuf::Descriptor* EntryList::descriptor() {
  protobuf_AssignDescriptorsOnce();
  return EntryList_descriptor_;
}

const EntryList& EntryList::default_instance() {
  if (default_instance_ == NULL) protobuf_AddDesc_index_2eproto();  return *default_instance_;
}

EntryList* EntryList::default_instance_ = NULL;

EntryList* EntryList::New() const {
  return new EntryList;
}

void EntryList::Clear() {
  entry_.Clear();
  ::memset(_has_bits_, 0, sizeof(_has_bits_));
  mutable_unknown_fields()->Clear();
}

bool EntryList::MergePartialFromCodedStream(
    ::google::protobuf::io::CodedInputStream* input) {
#define DO_(EXPRESSION) if (!(EXPRESSION)) return false
  ::google::protobuf::uint32 tag;
  while ((tag = input->ReadTag()) != 0) {
    switch (::google::protobuf::internal::WireFormatLite::GetTagFieldNumber(tag)) {
      // repeated .idxfile.Entry entry = 1;
      case 1: {
        if (::google::protobuf::internal::WireFormatLite::GetTagWireType(tag) ==
            ::google::protobuf::internal::WireFormatLite::WIRETYPE_LENGTH_DELIMITED) {
         parse_entry:
          DO_(::google::protobuf::internal::WireFormatLite::ReadMessageNoVirtual(
                input, add_entry()));
        } else {
          goto handle_uninterpreted;
        }
        if (input->ExpectTag(10)) goto parse_entry;
        if (input->ExpectAtEnd()) return true;
        break;
      }
      
      default: {
      handle_uninterpreted:
        if (::google::protobuf::internal::WireFormatLite::GetTagWireType(tag) ==
            ::google::protobuf::internal::WireFormatLite::WIRETYPE_END_GROUP) {
          return true;
        }
        DO_(::google::protobuf::internal::WireFormat::SkipField(
              input, tag, mutable_unknown_fields()));
        break;
      }
    }
  }
  return true;
#undef DO_
}

void EntryList::SerializeWithCachedSizes(
    ::google::protobuf::io::CodedOutputStream* output) const {
  // repeated .idxfile.Entry entry = 1;
  for (int i = 0; i < this->entry_size(); i++) {
    ::google::protobuf::internal::WireFormatLite::WriteMessageMaybeToArray(
      1, this->entry(i), output);
  }
  
  if (!unknown_fields().empty()) {
    ::google::protobuf::internal::WireFormat::SerializeUnknownFields(
        unknown_fields(), output);
  }
}

::google::protobuf::uint8* EntryList::SerializeWithCachedSizesToArray(
    ::google::protobuf::uint8* target) const {
  // repeated .idxfile.Entry entry = 1;
  for (int i = 0; i < this->entry_size(); i++) {
    target = ::google::protobuf::internal::WireFormatLite::
      WriteMessageNoVirtualToArray(
        1, this->entry(i), target);
  }
  
  if (!unknown_fields().empty()) {
    target = ::google::protobuf::internal::WireFormat::SerializeUnknownFieldsToArray(
        unknown_fields(), target);
  }
  return target;
}

int EntryList::ByteSize() const {
  int total_size = 0;
  
  // repeated .idxfile.Entry entry = 1;
  total_size += 1 * this->entry_size();
  for (int i = 0; i < this->entry_size(); i++) {
    total_size +=
      ::google::protobuf::internal::WireFormatLite::MessageSizeNoVirtual(
        this->entry(i));
  }
  
  if (!unknown_fields().empty()) {
    total_size +=
      ::google::protobuf::internal::WireFormat::ComputeUnknownFieldsSize(
        unknown_fields());
  }
  GOOGLE_SAFE_CONCURRENT_WRITES_BEGIN();
  _cached_size_ = total_size;
  GOOGLE_SAFE_CONCURRENT_WRITES_END();
  return total_size;
}

void EntryList::MergeFrom(const ::google::protobuf::Message& from) {
  GOOGLE_CHECK_NE(&from, this);
  const EntryList* source =
    ::google::protobuf::internal::dynamic_cast_if_available<const EntryList*>(
      &from);
  if (source == NULL) {
    ::google::protobuf::internal::ReflectionOps::Merge(from, this);
  } else {
    MergeFrom(*source);
  }
}

void EntryList::MergeFrom(const EntryList& from) {
  GOOGLE_CHECK_NE(&from, this);
  entry_.MergeFrom(from.entry_);
  mutable_unknown_fields()->MergeFrom(from.unknown_fields());
}

void EntryList::CopyFrom(const ::google::protobuf::Message& from) {
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

void EntryList::CopyFrom(const EntryList& from) {
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool EntryList::IsInitialized() const {
  
  for (int i = 0; i < entry_size(); i++) {
    if (!this->entry(i).IsInitialized()) return false;
  }
  return true;
}

void EntryList::Swap(EntryList* other) {
  if (other != this) {
    entry_.Swap(&other->entry_);
    std::swap(_has_bits_[0], other->_has_bits_[0]);
    _unknown_fields_.Swap(&other->_unknown_fields_);
    std::swap(_cached_size_, other->_cached_size_);
  }
}

::google::protobuf::Metadata EntryList::GetMetadata() const {
  protobuf_AssignDescriptorsOnce();
  ::google::protobuf::Metadata metadata;
  metadata.descriptor = EntryList_descriptor_;
  metadata.reflection = EntryList_reflection_;
  return metadata;
}


// @@protoc_insertion_point(namespace_scope)

}  // namespace idxfile

// @@protoc_insertion_point(global_scope)