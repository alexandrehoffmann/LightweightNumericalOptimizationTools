#include <gtest/gtest.h>

#include "LNOT/Containers/CircularBuffer.hpp"

// NOLINTBEGIN(cppcoreguidelines-avoid-magic-numbers) : Irrelevant
TEST(CircularBufferTests, basicOperations)
{
  lnot::CircularBuffer<int> buffer(3);
  EXPECT_EQ(buffer.getSize(), 0);
  buffer.push(1);
  EXPECT_EQ(buffer.getSize(), 1);
  buffer.push(2);
  EXPECT_EQ(buffer.getSize(), 2);
  buffer.push(3);
  EXPECT_EQ(buffer.getSize(), 3);
  buffer.push(4);
  EXPECT_EQ(buffer.getSize(), 3);
  EXPECT_EQ(*buffer.begin(), 2);
  EXPECT_EQ(*std::next(buffer.begin(), 1), 3);
  EXPECT_EQ(*std::next(buffer.begin(), 2), 4);
}

TEST(CircularBufferTests, iteratorOperations)
{
  lnot::CircularBuffer<int> buffer(3);
  buffer.push(1);
  buffer.push(2);
  buffer.push(3);
  auto it = buffer.begin();
  EXPECT_EQ(*it, 1);
  ++it;
  EXPECT_EQ(*it, 2);
  it++;
  EXPECT_EQ(*it, 3);
  it = buffer.begin();
  std::advance(it, 2);
  EXPECT_EQ(*it, 3);

  std::advance(it, -2);
  EXPECT_EQ(*it, 1);
}

TEST(CircularBufferTests, emptyBuffer)
{
  lnot::CircularBuffer<int> buffer(3);
  EXPECT_EQ(buffer.getSize(), 0);
  EXPECT_EQ(buffer.begin(), buffer.end());
}

TEST(CircularBufferTests, overwriteOldest)
{
  lnot::CircularBuffer<int> buffer(3);
  buffer.push(1);
  buffer.push(2);
  buffer.push(3);
  buffer.push(4);
  EXPECT_EQ(buffer.getSize(), 3);
  EXPECT_EQ(*buffer.begin(), 2);
  EXPECT_EQ(*std::next(buffer.begin(), 1), 3);
  EXPECT_EQ(*std::next(buffer.begin(), 2), 4);

  const auto bound = std::ranges::next(buffer.begin(), buffer.end());
  EXPECT_EQ(bound, buffer.end());
  EXPECT_EQ(std::prev(bound), std::next(buffer.begin(), 2));
  EXPECT_EQ(*std::prev(bound), 4);
  EXPECT_EQ(std::prev(bound, 2), std::next(buffer.begin(), 1));
  EXPECT_EQ(*std::prev(bound, 2), 3);
  EXPECT_EQ(std::prev(bound, 3), buffer.begin());
  EXPECT_EQ(*std::prev(bound, 3), 2);
}

TEST(CircularBufferTests, pushAndIterate)
{
  lnot::CircularBuffer<int> buffer(3);
  buffer.push(1);
  buffer.push(2);
  buffer.push(3);
  buffer.push(4);
  int expected = 2;
  for (const int& value : buffer)
  {
    EXPECT_EQ(value, expected);
    ++expected;
  }
}

TEST(CircularBufferTests, clearBuffer)
{
  lnot::CircularBuffer<int> buffer(3);
  buffer.push(1);
  buffer.push(2);
  buffer.push(3);
  buffer.clear();
  EXPECT_EQ(buffer.getSize(), 0);
  EXPECT_EQ(buffer.begin(), buffer.end());
}

TEST(CircularBufferTests, pushAfterClear)
{
  lnot::CircularBuffer<int> buffer(3);
  buffer.push(1);
  buffer.push(2);
  buffer.push(3);
  buffer.clear();
  buffer.push(4);
  EXPECT_EQ(buffer.getSize(), 1);
  EXPECT_EQ(*buffer.begin(), 4);
}

TEST(CircularBufferTests, pushAfterClearAndOverwrite)
{
  lnot::CircularBuffer<int> buffer(3);
  buffer.push(1);
  buffer.push(2);
  buffer.push(3);
  buffer.clear();
  buffer.push(4);
  buffer.push(5);
  buffer.push(6);
  buffer.push(7);
  EXPECT_EQ(buffer.getSize(), 3);
  EXPECT_EQ(*buffer.begin(), 5);
  EXPECT_EQ(*std::next(buffer.begin(), 1), 6);
  EXPECT_EQ(*std::next(buffer.begin(), 2), 7);
}

TEST(CircularBufferTests, offsetAfterOverwrite)
{
  lnot::CircularBuffer<int> buffer(3);
  buffer.push(1);
  buffer.push(2);
  buffer.push(3);
  buffer.push(4);
  EXPECT_EQ(buffer.getBoundOffset(), 1);
}

TEST(CircularBufferTests, offsetAfterClear)
{
  lnot::CircularBuffer<int> buffer(3);
  buffer.push(1);
  buffer.push(2);
  buffer.push(3);
  buffer.clear();
  EXPECT_EQ(buffer.getBoundOffset(), 0);
}

TEST(CircularBufferTests, offsetAfterPush)
{
  lnot::CircularBuffer<int> buffer(3);
  buffer.push(1);
  EXPECT_EQ(buffer.getBoundOffset(), 1);
  buffer.push(2);
  EXPECT_EQ(buffer.getBoundOffset(), 2);
  buffer.push(3);
  EXPECT_EQ(buffer.getBoundOffset(), 0);
  buffer.push(4);
  EXPECT_EQ(buffer.getBoundOffset(), 1);
}

TEST(CircularBufferTests, offsetAfterClearAndPush)
{
  lnot::CircularBuffer<int> buffer(3);
  buffer.push(1);
  buffer.push(2);
  buffer.push(3);
  buffer.clear();
  EXPECT_EQ(buffer.getBoundOffset(), 0);
  buffer.push(4);
  EXPECT_EQ(buffer.getBoundOffset(), 1);
  buffer.push(5);
  EXPECT_EQ(buffer.getBoundOffset(), 2);
  buffer.push(6);
  EXPECT_EQ(buffer.getBoundOffset(), 0);
}
// NOLINTEND(cppcoreguidelines-avoid-magic-numbers)
