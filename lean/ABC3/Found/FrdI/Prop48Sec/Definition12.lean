/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop44Gl
import ABC3.Found.FrdI.Def27
import ABC3.Found.FrdI.Prop48Sec.Proposition48

/-!
# Prop48Sec —— `[FrdI] Definition 1.2` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.FrdI

universe v u w u2 v2
open CategoryTheory
variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} (G : Frobenioid P)

/-! ## ★★★`iii-b-compact` の第 2 条 —— 無限位数は**逆向き**に落ちる

原文 (FrdI p.23):
> C such that for every n ∈N≥1, it holds that every B ∈Ob(C) base-isomorphic to

★★`Frobenius-compact` の 3 条のうち第 2 条は
「**無限位数の単元が存在する**」である。

★★★**測って分かったこと(2026-08-19)**: この条だけは
**準同型の逆向きに落ちる** —— `f : M →* N` と `u : N` が無限位数、
`v : M` が `f v = u` を満たすなら、`v` も無限位数。
★`v^k = 1` なら `f (v^k) = u^k = 1` で矛盾するだけ。

★★これが `iii-b-compact` の隔たり(`((𝒞^istr)^un-tr)^birat` の compact 対象を
`(𝒞^istr)^birat` へ渡す)に効く ——
**`𝒪^×` の写像が全射でありさえすれば第 2 条は渡る**。
★残る第 1 条(可換性)と第 3 条(`c/d` 倍作用)は向きが違うので、別に測る。 -/

/-- ★★★★**無限位数は準同型の逆向きに落ちる**。

★`v^k = 1` なら像も `1` になるので、像が無限位数なら `v` も無限位数。 -/
theorem infiniteOrder_of_map {M N : Type*} [Monoid M] [Monoid N] (f : M →* N)
    {u : N} (hu : ∀ k : ℕ+, u ^ ((k : ℕ+) : ℕ) ≠ 1)
    {v : M} (hv : f v = u) (k : ℕ+) : v ^ ((k : ℕ+) : ℕ) ≠ 1 := by
  intro hk
  refine hu k ?_
  rw [← hv, ← map_pow, hk, map_one]

def infiniteOrder_of_map.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 23,
    item := "Definition 1.2, (iv) — Frobenius-compact の第 2 条は逆向きに落ちる",
    sectionId := "frdi-def-1-2-iv" }

end ABC3.Found.FrdI
