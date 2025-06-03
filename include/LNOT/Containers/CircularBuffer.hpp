#ifndef LNOT_CIRCULAR_BUFFER_HPP
#define LNOT_CIRCULAR_BUFFER_HPP

#include <memory>

/**
 * @class CircularBuffer
 * 
 * Implements a fixed size fifo buffer.
 * Elements are stored in a standard C array.
 * Elements are added via the push method that returns the index at wich the elements is stored.
 * The buffer can be traversed in order via the foreach / reverseForeach methods : 
 * @code
 * buffer.foreach([](const size_type storage_index, const value_type& element_value)
 * {
 * 		// do stuff
 * });
 * @endcode
 * or with standard C++ iterators
 * @code
 * for (const value_type& element_value : buffer)
 * {
 * 		// do stuff
 * }
 * 
 * for (const_iterator it=buffer.begin(); it!=buffer.end(); ++it)
 * {
 * 		const size_type   storage_index = buffer.getIndex(it);
 *    const value_type& element_value = *it;
 * }
 * @endcode
 * 
 * This buffer is used by the LBFGS solver.
 * 
 */

namespace LNOT
{
	
template<typename T, typename Allocator = std::allocator<T> >
class CircularBuffer
{
	using AllocatorTraits = std::allocator_traits<Allocator>;

	template<bool isConst> class Iterator;
	template<bool isConst> class ReverseIterator;

	struct         sentinel {};
	struct reverse_sentinel {};
public:
	using value_type             = typename AllocatorTraits::value_type;
	using allocator              = typename AllocatorTraits::allocator_type;
	using size_type              = typename AllocatorTraits::size_type;
	using difference_type        = typename AllocatorTraits::difference_type;
	using reference              = value_type &;
	using const_reference        = const value_type &;
	using pointer                = typename AllocatorTraits::pointer;
	using const_pointer          = typename AllocatorTraits::const_pointer;
	
	static_assert (std::is_same<T, value_type>::value);
private:
	template<typename Function> struct IsConstTraverserFunction : std::bool_constant<std::is_invocable<Function, size_type, const_reference>::value> {};
	template<typename Function> struct IsTraverserFunction      : std::bool_constant<std::is_invocable<Function, size_type,       reference>::value> {};
public:
	CircularBuffer() = default;
	CircularBuffer(const size_type size, const allocator& alloc = allocator());
	CircularBuffer(const CircularBuffer &other);
	
	~CircularBuffer();

	size_type size() const { return size_type(m_data_end - m_data_begin); }

	bool allFilled() const { return m_allFilled; }

	template<typename... Args> 
	size_type push(Args &&...args);
	size_type push(const_reference value) { return push(value); }
	
	size_type firstInsertedIndex() const { return (m_allFilled and m_range_end!=m_data_end) ? size_type(m_range_end - m_data_begin) : 0u; }
	
	template<typename Function> Function foreach(Function f)              requires (IsTraverserFunction<Function>::value);
	template<typename Function> Function reverseForeach(Function f)       requires (IsTraverserFunction<Function>::value);
	
	template<typename Function> Function foreach(Function f)        const requires (IsConstTraverserFunction<Function>::value);
	template<typename Function> Function reverseForeach(Function f) const requires (IsConstTraverserFunction<Function>::value);

	void clear() { m_range_end = m_data_begin; m_allFilled = false; }
	
	bool empty() const { return (not m_allFilled) and m_range_end == m_data_begin; }
private:	
	allocator m_allocator;
	pointer   m_data_begin = nullptr;
	pointer   m_data_end   = nullptr;
	pointer   m_range_end  = nullptr;
	bool      m_allFilled  = false;
};
	
} // namespace LightOptim

#include <LNOT/Containers/CircularBuffer_impl.hpp>

#endif // LNOT_CIRCULAR_BUFFER_HPP
