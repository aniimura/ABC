import ABC3.Found.FrdI.Prop44Pre
import ABC3.Found.FrdI.Def27
import ABC3.Found.FrdI.Prop113

/-!
# [FrdI] Definition 4.5 —— 実装できる 2 条(と、待っている 2 条)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.86。

原文 (FrdI p.86):
> Definition 4.5.

## ★★4 条の依存(測定、2026-08-17)

| 条 | 内容 | 依存 | 状況 |
|---|---|---|---|
| (i) | `birationally Frobenius-normalized` / 同型 / `model type` | ★**`𝒞^birat` の pre-Frobenioid 構造**(取得済み) | ★**本ファイルで実装** |
| (ii) | `strictly rational` / `rational` / 同型 | ★★**`Φ^birat`**(`Proposition 4.4, (iii)`) | 待ち |
| (iii) | `rationally standard type` | (ii) ＋ `(𝒞^un-tr)^birat` ＋ `Frobenius-compact` | 待ち |
| (iv) | `Div-slim` | ★**`Φ` の関手性だけ** | ★**本ファイルで実装** |

★★**(iv) はこの命題群の中で唯一「`𝒞` を見ない」条件**である ——
`𝒟` と `Φ` だけで決まる。★`Corollary 4.11` / `Corollary 5.4` / `Corollary 5.7` が
いずれも仮定として引くので、**先に置いておく価値がある**。

## ★(iv) の Lean での書き方

原文 (FrdI p.86):
> Aut(DA →D) →Aut(DA →Mon)

★我々の `MonoidOn` は `Φ.functor : 𝒟ᵒᵖ ⥤ AddCommMon` という**反変**の形なので、
`𝒟_A → Mon` は **`(𝒟_A)ᵒᵖ ⥤ AddCommMon`**、すなわち
`(Over.forget A).op ⋙ Φ.functor` である。
★`Aut` の間の写像は `NatIso.op` と `isoWhiskerRight` の合成で書ける。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D}

/-! ## ★(iv) —— `Div-slim` -/

variable (Φ) in
/-- ★**`Φ` をスライス `𝒟_A` の上へ引いた関手** —— 原文の `𝒟_A → Mon`。

★`Φ` は反変なので、行き先は `(𝒟_A)ᵒᵖ ⥤ AddCommMon` である。 -/
def overPhi (A : D) : (Over A)ᵒᵖ ⥤ AddCommMonCat.{w} :=
  (Over.forget A).op ⋙ Φ.functor

variable (Φ) in
/-- ★★**`Aut(𝒟_A → 𝒟) → Aut(𝒟_A → Mon)`** —— 原文の「induced by composition with
the functor `Φ`」。 -/
def overPhiAut (A : D) (η : Aut (Over.forget A)) : Aut (overPhi Φ A) :=
  Functor.isoWhiskerRight (NatIso.op η) Φ.functor

variable (Φ) in
/-- ★★★**[FrdI] Definition 4.5, (iv)** —— `𝒟` が `Φ` に関して **Div-slim**。

原文 (FrdI p.86):
> Aut(DA →D) →Aut(DA →Mon)

★★**`𝒞` を一切見ない条件**である —— `𝒟` と `Φ` だけで決まる。 -/
def IsDivSlim : Prop := ∀ A : D, Function.Injective (overPhiAut Φ A)

variable (Φ) in
/-- ★★**slim なら Div-slim**(原文の「Thus, if `D` is slim, then it is Div-slim.」)。

原文 (FrdI p.86):
> Div-slim [relative to Φ] if, for every A ∈Ob(D),

★`Aut(𝒟_A → 𝒟)` が自明なら、そこからの写像は何であれ単射である。 -/
theorem isDivSlim_of_isSlim (hslim : IsSlimCat D) : IsDivSlim Φ := by
  intro A η₁ η₂ _
  rw [hslim A η₁, hslim A η₂]

/-! ## ★(i) —— `birationally Frobenius-normalized` -/

variable (P : PreFrobenioid C Φ) (G : Frobenioid P)

/-- ★★★**[FrdI] Definition 4.5, (i)** —— `𝒞` の対象が
**birationally Frobenius-normalized** であること。

原文 (FrdI p.86):
> (i) We shall say that an object of C is birationally Frobenius-normalized if

★★`𝒞^birat` の pre-Frobenioid 構造(`biratPre`)の上で
`Frobenius-normalized` であることをそのまま言う。 -/
def IsBirationallyFrobeniusNormalized (A : C) : Prop :=
  IsFrobeniusNormalized (biratPre P G) ((toBiratCat P G).obj A)

variable (C) in
/-- ★★**[FrdI] Definition 4.5, (i)** —— `𝒞` が
**birationally Frobenius-normalized 型**であること。 -/
def IsOfBirationallyFrobeniusNormalizedType : Prop :=
  ∀ A : C, IsBirationallyFrobeniusNormalized P G A

variable (C) in
/-- ★★★**[FrdI] Definition 4.5, (i)** —— `𝒞` が **model 型**であること。

原文 (FrdI p.86):
> normalized type. If C is of pre-model and birationally Frobenius-

★`pre-model 型`(`Definition 2.7, (iii)`)かつ
`birationally Frobenius-normalized 型`。 -/
def IsOfModelType [IsConnected D] : Prop :=
  IsPreModelType P ∧ IsOfBirationallyFrobeniusNormalizedType C P G

/-! ## ★待っている 2 条 —— (ii) と (iii)

原文 (FrdI p.86):
> (ii) Suppose that Φ is perf-factorial; A ∈Ob(C). Then we shall say that A

★★(ii) の `strictly rational` は **`Φ^birat(A) ⊆ Φ^gp(A)`** の元
`a - b`(`a, b ∈ Φ(A)`、`p ∈ Supp(a)`、`p ∉ Supp(b)`)の存在を言う。
★`Φ^birat` は `Proposition 4.4, (iii)` が与えるもので、**未実装**である。

★★**測定**: `Prop44.lean` は `𝒞^birat → 𝔽_{0_𝒟}` までを作っており、
原文の図式の中段 **`𝒞^birat → 𝔽_{Φ^gp}`** は無い。
★`Φ^birat` はその中段の像として定まるので、**中段を作るのが先**である。
★中段の `Div` は代表元 `(a : A′ → A, φ : A′ → B)` に対し
`(Base a)^{-1}(Div φ) − (Base a)^{-1}(Div a) ∈ Φ^gp(A)` で与えられる(紙の上)。

★(iii) は (ii) に加えて `(𝒞^un-tr)^birat` と `Frobenius-compact` 対象が要る。
-/

end ABC3.Found.FrdI
