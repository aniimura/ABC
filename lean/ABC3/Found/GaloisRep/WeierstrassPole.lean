import ABC3.Found.GaloisRep.EllipticConst

/-!
# Galois (G6) 第 243 ブロック —— **★★★★格子点での主要部**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★見積もりを縮められた

§9-555 では葉 (c) の段 2(`℘` の加法定理)を **30–60 ブロック**と見積もった。
「零点の和 = 極の和」(Abel、偏角の原理)が要る、という読みだった。

★★しかし**Liouville だけで済む古典的な証明がある**:

    F(z) := ℘(z+w) − [ (1/4)·((℘'(z)−℘'(w))/(℘(z)−℘(w)))² − ℘(z) − ℘(w) ]

は `z` の楕円関数で、**見かけの極がすべて相殺する**。よって第 241・242 で定数、
値を一点で見れば 0 ——これが加法定理である。偏角の原理は要らない。

★★★極の相殺を確かめるには**主要部が要る**。それは mathlib に在った:

| 補題 | 内容 |
|---|---|
| `weierstrassPExcept_def` | `℘[L−l₀](z) = ℘(z) + (1/l₀² − 1/(z−l₀)²)` |
| `analyticAt_weierstrassPExcept` | `℘[L−l₀]` は `l₀` で解析的 |
| `derivWeierstrassPExcept_def` | `℘'[L−l₀](z) = ℘'(z) + 2/(z−l₀)³` |
| `analyticAt_derivWeierstrassPExcept` | 同上 |

★つまり **`℘` から `l₀` の項だけ抜いたもの**が既に用意されていて、
それが `l₀` で解析的であることまで示されている。これで

    ℘(z)  = 1/(z−l₀)² + (正則部)
    ℘'(z) = −2/(z−l₀)³ + (正則部)

が直ちに書ける。★★見積もりは **15–35 ブロック**に縮まった。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `weierstrassP_principal` | ★★★★`℘` の主要部 |
| `derivWeierstrassP_principal` | ★★★★`℘'` の主要部 |
| `sq_mul_weierstrassP` | ★★★`(z−l₀)²℘(z) = 1 + O((z−l₀)²)` |
| `cube_mul_derivWeierstrassP` | ★★★`(z−l₀)³℘'(z) = −2 + O((z−l₀)³)` |
-/

namespace ABC3.Found.GaloisRep

open Complex Real PeriodPair

/-- ★★★★**格子点での `℘` の主要部** —— `℘(z) = 1/(z−l₀)² + (正則部)`。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem weierstrassP_principal (L : PeriodPair) (l₀ : L.lattice) (z : ℂ) :
    L.weierstrassP z
      = 1 / (z - (l₀ : ℂ)) ^ 2 + (L.weierstrassPExcept (l₀ : ℂ) z - 1 / (l₀ : ℂ) ^ 2) := by
  rw [L.weierstrassPExcept_def l₀ z]
  ring

/-- ★★★★**格子点での `℘'` の主要部** —— `℘'(z) = −2/(z−l₀)³ + (正則部)`。 -/
theorem derivWeierstrassP_principal (L : PeriodPair) (l₀ : L.lattice) (z : ℂ) :
    L.derivWeierstrassP z
      = -2 / (z - (l₀ : ℂ)) ^ 3 + L.derivWeierstrassPExcept (l₀ : ℂ) z := by
  rw [L.derivWeierstrassPExcept_def l₀ z]
  ring

/-- ★★★`(z−l₀)²·℘(z) = 1 + (z−l₀)²·(正則部)`。 -/
theorem sq_mul_weierstrassP (L : PeriodPair) (l₀ : L.lattice) (z : ℂ) (hz : z ≠ (l₀ : ℂ)) :
    (z - (l₀ : ℂ)) ^ 2 * L.weierstrassP z
      = 1 + (z - (l₀ : ℂ)) ^ 2 * (L.weierstrassPExcept (l₀ : ℂ) z - 1 / (l₀ : ℂ) ^ 2) := by
  have hne : (z - (l₀ : ℂ)) ≠ 0 := sub_ne_zero.2 hz
  rw [weierstrassP_principal L l₀ z]
  field_simp

/-- ★★★`(z−l₀)³·℘'(z) = −2 + (z−l₀)³·(正則部)`。 -/
theorem cube_mul_derivWeierstrassP (L : PeriodPair) (l₀ : L.lattice) (z : ℂ)
    (hz : z ≠ (l₀ : ℂ)) :
    (z - (l₀ : ℂ)) ^ 3 * L.derivWeierstrassP z
      = -2 + (z - (l₀ : ℂ)) ^ 3 * L.derivWeierstrassPExcept (l₀ : ℂ) z := by
  have hne : (z - (l₀ : ℂ)) ≠ 0 := sub_ne_zero.2 hz
  rw [derivWeierstrassP_principal L l₀ z]
  field_simp

/-! ## ★出典の紐付け(`.src`) -/

def weierstrassP_principal.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——格子点での ℘ の主要部)",
    sectionId := "genell-def-3-3" }

def derivWeierstrassP_principal.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——格子点での ℘' の主要部)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
