/*
This file incorporates code from the kPEX project.
The original source code is available at: https://github.com/ShuohaoGao/kPEX/blob/master/kPEX/MyBitset.h
*/
#ifndef MY_BITSET_H
#define MY_BITSET_H

#include "utility.h"

/**
 * @brief bitset, which supports & | bit-operation in O(n/64).
 * we provide an iterator too
 */
class MyBitset
{
public:
    int range;
    int n; // the size of buf is n+1
    int m; // valid size in last int (range & 63)
    uint64_t *buf;
    int sz;
    bool sz_changed; // lazy-tag

    MyBitset() : n(0), m(0), buf(nullptr), sz(0), sz_changed(false) {}

    MyBitset(int _range) : range(_range), n(_range >> 6), m(_range & 63), sz(0), sz_changed(false)
    {
        buf = new uint64_t[n + 1];
        memset(buf, 0, sizeof(uint64_t) * (n + 1));
    }

    MyBitset(const MyBitset &other) : range(other.range), n(other.n), m(other.m), sz(other.sz), sz_changed(other.sz_changed)
    {
        buf = new uint64_t[n + 1];
        memcpy(buf, other.buf, sizeof(uint64_t) * (n + 1));
    }

    ~MyBitset()
    {
        if (buf != nullptr) {
            delete[] buf;
            buf = nullptr;
        }
    }

    void clear()
    {
        memset(buf, 0, sizeof(uint64_t) * (n + 1));
        sz = 0;
        sz_changed = false;
    }

    void resize(int _range)
    {
        range = _range;
        n = _range >> 6;
        m = _range & 63;
        if (buf != nullptr) {
            delete[] buf;
            buf = nullptr;
        }
        buf = new uint64_t[n + 1];
        memset(buf, 0, sizeof(uint64_t) * (n + 1));
    }

    void flip()
    {
        for (int i = 0; i < n; ++i) buf[i] = ~buf[i];
        buf[n] ^= (1ULL << m) - 1;
        sz_changed = true;
    }

    inline void insert(int x)
    {
        assert(x < range);
        assert(!contains(x));
        buf[x >> 6] |= 1ULL << (x & 63);
        sz++;
    }

    inline void remove(int x)
    {
        assert(x < range);
        assert(contains(x));
        buf[x >> 6] &= ~(1ULL << (x & 63));
        sz--;
    }

    inline bool contains(int x) const
    {
        assert(x < range);
        return (buf[x >> 6] >> (x & 63)) & 1ULL;
    }

    bool empty() const
    {
        return std::all_of(buf, buf + n + 1, [](uint64_t val)
            { return val == 0ULL; });
    }

    MyBitset &operator=(const MyBitset &other)
    {
        if (this == &other) return *this;
        if (buf != nullptr) {
            delete[] buf;
            buf = nullptr;
        }
        range = other.range;
        n = other.n;
        m = other.m;
        sz = other.sz;
        sz_changed = other.sz_changed;
        buf = new uint64_t[n + 1];
        memcpy(buf, other.buf, sizeof(uint64_t) * (n + 1));
        return *this;
    }

    void operator&=(const MyBitset &other)
    {
        assert(n == other.n);
        std::transform(buf, buf + n + 1, other.buf, buf, std::bit_and<uint64_t>());
        sz_changed = true;
    }

    void operator|=(const MyBitset &other)
    {
        assert(n == other.n);
        std::transform(buf, buf + n + 1, other.buf, buf, std::bit_or<uint64_t>());
        sz_changed = true;
    }

    void operator^=(const MyBitset &other)
    {
        assert(n == other.n);
        std::transform(buf, buf + n + 1, other.buf, buf, std::bit_xor<uint64_t>());
        sz_changed = true;
    }

    bool operator==(const MyBitset &other) const
    {
        for (int i = 0; i <= n; i++) if (buf[i] != other.buf[i]) return false;
        return true;
    }

    int xor_size(const MyBitset &other) const
    {
        return std::inner_product(buf, buf + n + 1, other.buf, 0, std::plus<int>(), [](uint64_t a, uint64_t b)
            { return __builtin_popcountll(a ^ b); });
    }

    int intersect_size(const MyBitset &other) const
    {
        return std::inner_product(buf, buf + n + 1, other.buf, 0, std::plus<int>(), [](uint64_t a, uint64_t b)
            { return __builtin_popcountll(a & b); });
    }

    int union_size(const MyBitset &other) const
    {
        return std::inner_product(buf, buf + n + 1, other.buf, 0, std::plus<int>(), [](uint64_t a, uint64_t b)
            { return __builtin_popcountll(a | b); });
    }

    MyBitset operator|(const MyBitset &other) const
    {
        assert(range == other.range);
        MyBitset result(range);
        std::transform(buf, buf + n + 1, other.buf, result.buf, std::bit_or<uint64_t>());
        result.sz_changed = true;
        return result;
    }

    MyBitset operator&(const MyBitset &other) const
    {
        assert(range == other.range);
        MyBitset result(range);
        std::transform(buf, buf + n + 1, other.buf, result.buf, std::bit_and<uint64_t>());
        result.sz_changed = true;
        return result;
    }

    bool operator[](int x) const
    {
        assert(x < range);
        return (buf[x >> 6] >> (x & 63)) & 1ULL;
    }
    // a.sub(b) <==> a= a - a ∩ b
    void sub(const MyBitset &other)
    {
        for (int i = 0; i <= n; i++) buf[i] &= ~other.buf[i];
        sz_changed = true;
    }

    int size()
    {
        if (!sz_changed) return sz;
        sz = std::accumulate(buf, buf + n + 1, 0, [](int sum, uint64_t val)
            { return sum + __builtin_popcountll(val); });
        sz_changed = false;
        return sz;
    }

    int intersect(const MyBitset &other) const
    {
        return std::inner_product(buf, buf + n + 1, other.buf, 0, std::plus<int>(), [](uint64_t a, uint64_t b)
            { return __builtin_popcountll(a & b); });
    }

    int intersect(const MyBitset &a, const MyBitset &b) const
    {
        int ret = 0;
        for (int i = 0; i <= n; i++) ret += __builtin_popcountll((buf[i] & a.buf[i]) & b.buf[i]);
        return ret;
    }

    class Iterator
    {
    private:
        int index;
        const MyBitset &bitset;

    public:
        Iterator(int idx, const MyBitset &bs) : index(idx), bitset(bs) {}

        Iterator &operator++()
        {
            if (index < bitset.range - 1) {
                index++;
                int block_index = index >> 6;
                uint64_t mask = bitset.buf[block_index] >> (index & 63);
                while (mask == 0 && (block_index) < bitset.n) {
                    mask = bitset.buf[++block_index];
                    index = block_index << 6; // Move to the next block of 64 bits
                }
                if (mask != 0) {
                    index += __builtin_ctzll(mask); // Find the next set bit using __builtin_ctzll
                }
                else
                    index = bitset.range;
            }
            else {
                index = bitset.range;
            }
            return *this;
        }

        int operator*() const
        {
            return index;
        }

        bool operator!=(const Iterator &other) const
        {
            return index != other.index;
        }
    };

    Iterator begin() const
    {
        int idx = 0;
        int blockIndex = 0;
        uint64_t mask = buf[blockIndex];
        while (mask == 0 && blockIndex < n) {
            mask = buf[++blockIndex];
            idx += 64; // Move to the next block of 64 bits
        }
        if (mask != 0) {
            idx += __builtin_ctzll(mask); // Find the first set bit using __builtin_ctzll
        }
        else
            idx = range;
        return Iterator(idx, *this);
    }

    Iterator end() const
    {
        return Iterator(range, *this);
    }
};

#endif