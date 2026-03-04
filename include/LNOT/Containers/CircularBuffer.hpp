#ifndef LNOT_CIRCULAR_BUFFER_HPP
#define LNOT_CIRCULAR_BUFFER_HPP

#include <memory>

namespace LNOT
{

/**
 * @brief A fixed-size circular buffer container with customizable allocator.
 * 
 * This container manages a pre-allocated memory buffer and allows pushing elements in a circular fashion.
 * Once full, the buffer wraps around and overwrites the oldest values. Internally, it maintains a pointer
 * to the start and end of the data range, and supports iteration in forward and reverse directions.
 * 
 * Elements are inserted using `push()`, and traversal is supported via `foreach()` and `reverseForeach()`
 * using user-defined callable objects.
 * 
 * @tparam T Type of the elements stored in the buffer.
 * @tparam Allocator Allocator type used for memory management (defaults to `std::allocator<T>`).
 */	
template<typename T, typename Allocator = std::allocator<T> >
class CircularBuffer
{
	using AllocatorTraits = std::allocator_traits<Allocator>;
public:
	using value_type       = typename AllocatorTraits::value_type;      ///<  @brief Type of elements stored in the buffer.
	using allocator        = typename AllocatorTraits::allocator_type;  ///<  @brief Allocator type used to manage memory.
	using size_type        = typename AllocatorTraits::size_type;       ///<  @brief Type used to represent sizes.
	using difference_type  = typename AllocatorTraits::difference_type; ///<  @brief Type used to represent differences between pointers.
	using reference        = value_type &;                              ///<  @brief Reference to an element.
	using const_reference  = const value_type &;                        ///<  @brief Const reference to an element.
	using pointer          = typename AllocatorTraits::pointer;         ///<  @brief Pointer to an element.
	using const_pointer    = typename AllocatorTraits::const_pointer;   ///<  @brief Const pointer to an element.
	
	static_assert (std::is_same<T, value_type>::value);
private:
	template<typename Function> static constexpr bool isConstTraverserFunction = std::invocable<Function, size_type, const_reference>; 
	template<typename Function> static constexpr bool isTraverserFunction      = std::invocable<Function, size_type,       reference>;
public:
	CircularBuffer() = default;
	
	/**
	 * @brief Construct a circular buffer with a fixed size.
	 * 
	 * @param size Number of elements the buffer can hold.
	 * @param alloc Optional allocator instance.
	 */
	CircularBuffer(const size_type size, const allocator& alloc = allocator());
	
	/**
	 * @brief Copy constructor.
	 * 
	 * Performs a deep copy of the buffer content.
	 * 
	 * @param other Buffer to copy from.
	 */
	CircularBuffer(const CircularBuffer &other);
	
	~CircularBuffer();
	
	/// @brief Get the number of elements stored in the buffer.
	size_type size() const { return size_type(m_data_end - m_data_begin); }
	
	/// @brief Check if the buffer has been completely filled at least once.
	bool allFilled() const { return m_allFilled; }
	
	/**
	 * @brief Insert a new element in-place using constructor arguments.
	 * 
	 * @tparam Args Argument types for element construction.
	 * @param args Arguments forwarded to the constructor of T.
	 * @return The index (logical) of the newly inserted element.
	 */
	template<typename... Args> 
	size_type push(Args&&...args);
	
	/// @brief Get the index of the first inserted element relative to the buffer's internal memory.
	size_type firstInsertedIndex() const { return (m_allFilled and m_range_end!=m_data_end) ? size_type(m_range_end - m_data_begin) : 0u; }
	
	/**
	 * @brief Traverse all elements in order and apply a user function.
	 * 
	 * The function must accept (index, reference).
	 * 
	 * @tparam Function A callable that takes (size_type, reference).
	 * @param f The function to apply to each element.
	 * @return The same function object (may be used for stateful operations).
	 */
	template<typename Function> Function foreach(Function f) requires (isTraverserFunction<Function>);
	
	/**
	 * @brief Traverse all elements in reverse order and apply a user function.
	 * 
	 * The function must accept (index, reference).
	 * 
	 * @tparam Function A callable that takes (size_type, reference).
	 * @param f The function to apply to each element in reverse.
	 * @return The same function object.
	 */
	template<typename Function> Function reverseForeach(Function f) requires (isTraverserFunction<Function>);
	
	/**
	 * @brief Traverse all elements (const version) and apply a user function.
	 * 
	 * The function must accept (index, const_reference).
	 * 
	 * @tparam Function A callable that takes (size_type, const_reference).
	 * @param f The function to apply to each element.
	 * @return The same function object.
	 */
	template<typename Function> Function foreach(Function f) const requires (isConstTraverserFunction<Function>);
	
	/**
	 * @brief Traverse all elements in reverse (const version) and apply a user function.
	 * 
	 * The function must accept (index, const_reference).
	 * 
	 * @tparam Function A callable that takes (size_type, const_reference).
	 * @param f The function to apply in reverse order.
	 * @return The same function object.
	 */
	template<typename Function> Function reverseForeach(Function f) const requires (isConstTraverserFunction<Function>);
	
	/// @brief Clears the buffer contents (logically, does not deallocate).
	void clear() { m_range_end = m_data_begin; m_allFilled = false; }
	
	/// @brief Check if the buffer currently holds no elements.
	bool empty() const { return (not m_allFilled) and m_range_end == m_data_begin; }
private:	
	allocator m_allocator;
	pointer   m_data_begin = nullptr;
	pointer   m_data_end   = nullptr;
	pointer   m_range_end  = nullptr;
	bool      m_allFilled  = false;
};

using CircularBufferF = CircularBuffer<float>;
using CircularBufferD = CircularBuffer<double>;
	
} // namespace LightOptim

#include <LNOT/Containers/CircularBuffer_impl.hpp>

#endif // LNOT_CIRCULAR_BUFFER_HPP
