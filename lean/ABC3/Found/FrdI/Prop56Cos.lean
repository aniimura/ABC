/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop56Sec

/-!
# [FrdI] Proposition 5.6 の証明中の商単系 `E_M = M \ E`

原文 (FrdI p.106):
> structure on the set of cosets EM

★原文はここで、`E ⊆ End_𝒞(A)`(base-isomorphism のなす部分単系)と
部分群 `M ⊆ 𝒪^×(A) ⊆ E` に対し、

* 「任意の `m ∈ M`, `f ∈ E` に対し `f · m = m' · f` なる `m' ∈ M` がある」

という条件の下で、**左剰余類の集合 `E_M := M \ E = {M · f}` に単系構造が入り、
全射準同型 `E ↠ E_M` が定まる**ことを使う。

★★**このファイルはその一般論(純粋に単系の代数)を実装する**。
`M_p ⊆ 𝒪^×(A)` の構成(原文 p.106、pro-`p` の議論)は別の葉である
(依存グラフの鎖 `prop56` の節点 `p56-Mp` / `p56-limit`)。

## ★条件の読み

原文の `M` は「部分群」だが、`End A` は単系なので
**「可逆で逆元も `M` に入る部分単系」**として書いた(`IsShiftStable.inv`)。
★`𝒪^×(A)` 自身が群であること(`otimesGroup`)は在庫にある。
-/

namespace ABC3.Found.FrdI

universe uu

variable {G : Type uu} [Monoid G]

/-- ★★**[FrdI] Proposition 5.6 の証明中の設定**(原文 p.106)——
部分単系 `M ⊆ E` が「群であり、`E` の元で左右を入れ替えられる」こと。 -/
structure IsShiftStable (E M : Submonoid G) : Prop where
  /-- `M ⊆ E`。 -/
  le : M ≤ E
  /-- `M` の元は可逆で、逆元も `M` に入る。 -/
  inv : ∀ x ∈ M, ∃ y ∈ M, x * y = 1 ∧ y * x = 1
  /-- 原文の「for any `m ∈ M`, `f ∈ E`, there exists `m' ∈ M` such that `f · m = m' · f`」。 -/
  shift : ∀ m ∈ M, ∀ f ∈ E, ∃ m' ∈ M, f * m = m' * f

namespace IsShiftStable

variable {E M : Submonoid G} (h : IsShiftStable E M)

/-- ★左剰余類 `M · f` の同値関係。 -/
def setoid : Setoid E where
  r f g := ∃ m ∈ M, (g : G) = m * f
  iseqv := {
    refl := fun f => ⟨1, M.one_mem, (one_mul _).symm⟩
    symm := by
      rintro f g ⟨m, hm, hg⟩
      obtain ⟨y, hy, _, hyx⟩ := h.inv m hm
      exact ⟨y, hy, by rw [hg, ← mul_assoc, hyx, one_mul]⟩
    trans := by
      rintro f g k ⟨m, hm, hg⟩ ⟨m', hm', hk⟩
      exact ⟨m' * m, M.mul_mem hm' hm, by rw [hk, hg, mul_assoc]⟩ }

/-- ★★**商単系 `E_M = M \ E`**(原文 p.106 の `EM`)。 -/
def Cos : Type uu := Quotient h.setoid

/-- ★射影 `E ↠ E_M`。 -/
def cosMk (f : E) : h.Cos := Quotient.mk h.setoid f

theorem cosMk_eq_cosMk {f g : E} : h.cosMk f = h.cosMk g ↔ ∃ m ∈ M, (g : G) = m * f :=
  ⟨fun hh => Quotient.exact hh, fun hh => Quotient.sound hh⟩

theorem cosMk_surjective : Function.Surjective h.cosMk :=
  fun x => Quotient.inductionOn x fun f => ⟨f, rfl⟩

/-- ★★★**商の積が定まる** —— ★ここで `shift` を使う。

`(m₁ f) · (m₂ g) = m₁ (f m₂) g = m₁ m₂' f g ∈ M · (f g)`。 -/
noncomputable instance : Monoid h.Cos where
  mul := Quotient.map₂ (fun f g : E => f * g) (by
    rintro f f' ⟨m, hm, hf⟩ g g' ⟨m', hm', hg⟩
    obtain ⟨m'', hm'', hshift⟩ := h.shift m' hm' f f.2
    refine ⟨m * m'', M.mul_mem hm hm'', ?_⟩
    show ((f' : G) * g') = m * m'' * (f * g)
    rw [hf, hg, mul_assoc, ← mul_assoc (f : G) m' (g : G), hshift]
    simp only [mul_assoc])
  one := h.cosMk 1
  one_mul x := Quotient.inductionOn x fun f => congrArg h.cosMk (one_mul f)
  mul_one x := Quotient.inductionOn x fun f => congrArg h.cosMk (mul_one f)
  mul_assoc x y z :=
    Quotient.inductionOn₃ x y z fun a b c => congrArg h.cosMk (mul_assoc a b c)

@[simp] theorem cosMk_mul (f g : E) : h.cosMk (f * g) = h.cosMk f * h.cosMk g := rfl

@[simp] theorem cosMk_one : h.cosMk 1 = 1 := rfl

/-- ★★**射影は単系の全射準同型** —— 原文の
「a natural monoid structure on the set of cosets `E_M`, together with a natural
surjection of monoids `E ↠ E_M`」。 -/
noncomputable def cosHom : E →* h.Cos where
  toFun := h.cosMk
  map_one' := rfl
  map_mul' _ _ := rfl

theorem cosHom_surjective : Function.Surjective h.cosHom := h.cosMk_surjective

/-- ★`M` の元は商で `1` に落ちる。 -/
theorem cosMk_eq_one_of_mem {f : E} (hf : (f : G) ∈ M) : h.cosMk f = 1 :=
  ((h.cosMk_eq_cosMk).mpr ⟨(f : G), hf, by simp⟩).symm

end IsShiftStable

/-! ### ★出典の紐付け -/

/-- ★locator —— `Proposition 5.6` の証明中の商単系 `E_M`(★**条つき**:
`M_p` の構成そのものは別の葉)。 -/
def IsShiftStable.cosHom.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 106,
    item := "Proposition 5.6 — 証明中の商単系 E_M = M \\ E と全射 E ↠ E_M",
    sectionId := "frdi-prop-5-6" }

end ABC3.Found.FrdI
