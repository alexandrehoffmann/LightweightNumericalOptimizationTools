#include <gtest/gtest.h>

#include "LNOT/Containers/ArrayHandle.hpp"

using ArrayHandle = lnot::HostArrayHandle<int>;

constexpr int SMALL_SIZE = 5;
constexpr int SIZE = 2 * SMALL_SIZE;

// NOLINTBEGIN(cppcoreguidelines-pro-bounds-avoid-unchecked-container-access) : Irrelevant
TEST(ArrayHandleTests, basicFunctionality)
{
  ArrayHandle arr = ArrayHandle::allocate(SMALL_SIZE);
  for (int i = 0; i != SMALL_SIZE; ++i)
  {
    arr[i] = i;
  }
  for (int i = 0; i != SMALL_SIZE; ++i)
  {
    EXPECT_EQ(arr[i], i);
  }
}

TEST(ArrayHandleTests, moveAssignment)
{
  ArrayHandle arr1 = ArrayHandle::allocate(SMALL_SIZE);
  for (int i = 0; i != SMALL_SIZE; ++i)
  {
    arr1[i] = i;
  }
  ArrayHandle arr2 = std::move(arr1);
  for (int i = 0; i != SMALL_SIZE; ++i)
  {
    EXPECT_EQ(arr2[i], i);
  }
  // NOLINTNEXTLINE(bugprone-use-after-move) : point of the test
  EXPECT_FALSE(static_cast<bool>(arr1));
}

TEST(ArrayHandleTests, moveAssignmentOperator)
{
  ArrayHandle arr1 = ArrayHandle::allocate(SMALL_SIZE);
  for (int i = 0; i != SMALL_SIZE; ++i)
  {
    arr1[i] = i;
  }
  ArrayHandle arr2 = std::move(arr1);
  for (int i = 0; i != SMALL_SIZE; ++i)
  {
    EXPECT_EQ(arr2[i], i);
  }
  // NOLINTNEXTLINE(bugprone-use-after-move) : point of the test
  EXPECT_FALSE(static_cast<bool>(arr1));
}

TEST(ArrayHandleTests, free)
{
  ArrayHandle arr = ArrayHandle::allocate(SMALL_SIZE);
  for (int i = 0; i != SMALL_SIZE; ++i)
  {
    arr[i] = i;
  }
  arr.free();
  EXPECT_FALSE(static_cast<bool>(arr));
}

TEST(ArrayHandleTests, operatorBool)
{
  ArrayHandle arr = ArrayHandle::allocate(SMALL_SIZE);
  EXPECT_TRUE(static_cast<bool>(arr));
  arr.free();
  EXPECT_FALSE(static_cast<bool>(arr));
}

TEST(ArrayHandleTests, operatorStar)
{
  ArrayHandle arr = ArrayHandle::allocate(SMALL_SIZE);
  for (int i = 0; i != SMALL_SIZE; ++i)
  {
    arr[i] = i;
  }
  EXPECT_EQ(*arr, 0);
}

TEST(ArrayHandleTests, reallocate)
{
  ArrayHandle arr = ArrayHandle::allocate(SMALL_SIZE);
  for (int i = 0; i != SMALL_SIZE; ++i)
  {
    arr[i] = i;
  }
  arr = ArrayHandle::allocate(SIZE);
  for (int i = 0; i != SIZE; ++i)
  {
    arr[i] = i;
  }
  for (int i = 0; i != SIZE; ++i)
  {
    EXPECT_EQ(arr[i], i);
  }
}

TEST(ArrayHandleTests, swap)
{
  ArrayHandle arr1 = ArrayHandle::allocate(SMALL_SIZE);
  ArrayHandle arr2 = ArrayHandle::allocate(SMALL_SIZE);
  for (int i = 0; i != SMALL_SIZE; ++i)
  {
    arr1[i] = i;
    arr2[i] = i + SMALL_SIZE;
  }
  arr1.swap(arr2);
  for (int i = 0; i != SMALL_SIZE; ++i)
  {
    EXPECT_EQ(arr1[i], i + SMALL_SIZE);
    EXPECT_EQ(arr2[i], i);
  }
}
// NOLINTEND(cppcoreguidelines-pro-bounds-avoid-unchecked-container-access)
