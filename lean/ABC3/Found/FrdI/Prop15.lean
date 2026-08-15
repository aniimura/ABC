import ABC3.Found.FrdI.Prop14

/-!
# [FrdI] Proposition 1.5 —— Elementary Frobenioids are Frobenioids

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、
物理 p.27(**400 dpi 目視確認 2026-08-15**)。

原文 (FrdI p.27):
> Proposition 1.5. (Elementary Frobenioids are Frobenioids) Let Φ be a

原文 (FrdI p.27):
> pre-divisorial monoid on a connected, totally epimorphic category D. Then:

原文 (FrdI p.27):
> (i) FΦ, equipped with the natural functor FΦ →FΦchar, is a Frobenioid of Aut-

## ★★`wP` との差 —— 何が模型固有だったか

`WitnessFrobenioid.lean` の `wIsFrobenioid` は `Definition 1.3` の22条を
`wP`(`𝒟 = Vee`、`Φ = ℕ`、`𝒞 = 𝔽_Φ`、`toElem = 𝟭`)について埋めた。
★**一般化にあたって、そのうち何が模型固有だったかを先に洗い出す。**

**模型固有だったもの(一般には使えない)**:

1. ★★**`wHom_ext`** —— 「`Vee` の hom は subsingleton なので射は `(wd, wn)` で決まる」。
   `wIsFrobenioid` の**ほぼ全条**がこれを使っている。一般の `𝒟` では
   底の射の等式を**別途示す**必要がある。
2. ★**`toElem = 𝟭`** —— `wP` は恒等関手だったので `Div φ = φ.div` だった。
   ★**原文の構造は `𝔽_Φ → 𝔽_{Φ^char}`** なので `Div φ = toChar φ.div` であり、
   **isometric は「零因子が 0」ではなく「零因子が可逆」を意味する**
   (`toChar_eq_zero_iff`)。`wP` では `Φ = ℕ` が sharp だったので両者が一致していた。

**模型固有でなかったもの(一般に成り立つ)**:

3. ★**`𝔽_Φ` は isotropic type** —— isometric な pre-step は
   「零因子が可逆・次数 1・底が同型」なので、**逆射を明示的に作れる**。
   下の `elemFrob_isotropic` で一般に証明する。
4. ★したがって **`𝔽_Φ` のすべての射が co-angular** ——
   これは `Proposition 1.4, (i)` の直接の帰結であり、
   ★**原文自身が同じ道を通る**:

原文 (FrdI p.27):
> sition 1.4, (i) [which is applicable to all morphisms of FΦ since FΦ is of isotropic

## ★このファイルの到達点

* `elemFrobToChar` —— 関手 `𝔽_Φ → 𝔽_{Φ^char}`
* `elemPreFrobenioid` —— ★**`𝔽_Φ` が `Φ^char` 上の pre-Frobenioid であること**
* `elemFrob_isotropic` / `elemFrob_coAngular` —— 上の 3, 4

★**`Definition 1.3` の22条は、まだ埋めていない。** どこで止まっているかは
ファイル末尾の記録を見よ。
-/

namespace ABC3.Found.FrdI

open CategoryTheory Opposite

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] (Φ : MonoidOn.{v, u, w} D)

/-! ### ★関手 `𝔽_Φ → 𝔽_{Φ^char}`

原文 (FrdI p.27):
> (i) FΦ, equipped with the natural functor FΦ →FΦchar, is a Frobenioid of Aut-
-/

/-- ★**`𝔽_Φ → 𝔽_{Φ^char}`** —— 対象はそのまま、零因子を `Φ(A)^char` へ送る。

関手性は `charMap_toChar`(`charMap g ∘ toChar = toChar ∘ g`)から出る。 -/
def elemFrobToChar : ElemFrobCat Φ ⥤ ElemFrobCat Φ.charOn where
  obj A := ⟨A.base⟩
  map {A B} f := ⟨f.base, toChar f.div, f.deg⟩
  map_id A := by
    refine ElemFrobCat.Hom.ext rfl ?_ rfl
    show toChar (0 : Φ.val A.base) = 0
    exact map_zero _
  map_comp {A B E} f g := by
    refine ElemFrobCat.Hom.ext rfl ?_ rfl
    show toChar (Φ.map f.base g.div + ((g.deg : ℕ+) : ℕ) • f.div)
      = charMap (Φ.map f.base) (toChar g.div) + ((g.deg : ℕ+) : ℕ) • toChar f.div
    rw [map_add, map_nsmul, charMap_toChar]

@[simp] theorem elemFrobToChar_div {A B : ElemFrobCat Φ} (f : A ⟶ B) :
    ((elemFrobToChar Φ).map f).div = toChar f.div := rfl
@[simp] theorem elemFrobToChar_deg {A B : ElemFrobCat Φ} (f : A ⟶ B) :
    ((elemFrobToChar Φ).map f).deg = f.deg := rfl
@[simp] theorem elemFrobToChar_base {A B : ElemFrobCat Φ} (f : A ⟶ B) :
    ((elemFrobToChar Φ).map f).base = f.base := rfl

/-! ### ★`𝔽_Φ` は `Φ^char` 上の pre-Frobenioid

原文 (FrdI p.27):
> (i)]; and the injectivity condition of Definition 1.1, (ii), (a). Thus, FΦ is a pre-
-/

variable (hD : IsTotallyEpimorphic D) (hpd : ∀ A : D, IsPreDivisorial (Φ.val A))

/-- ★**`Proposition 1.5, (i)` の前半** —— `𝔽_Φ` は `Φ^char` 上の pre-Frobenioid。

原文が挙げる3つの入力がそのままフィールドを埋める:

* `divisorial` ← `charOn_isDivisorialOn`(`Φ` が pre-divisorial なら `Φ^char` は divisorial)
* `totEpiC` ← `isTotallyEpimorphic_elemFrobCat`(`𝒟` が totally epimorphic + `Φ` が integral)
* `totEpiD` ← 仮定 -/
def elemPreFrobenioid : PreFrobenioid (ElemFrobCat Φ) Φ.charOn where
  toElem := elemFrobToChar Φ
  divisorial := Φ.charOn_isDivisorialOn hpd
  totEpiC := isTotallyEpimorphic_elemFrobCat hD (fun A => (hpd A).1)
  totEpiD := hD

@[simp] theorem elemPreFrobenioid_Div {A B : ElemFrobCat Φ} (f : A ⟶ B) :
    (elemPreFrobenioid Φ hD hpd).Div f = toChar f.div := rfl
@[simp] theorem elemPreFrobenioid_degFr {A B : ElemFrobCat Φ} (f : A ⟶ B) :
    (elemPreFrobenioid Φ hD hpd).degFr f = f.deg := rfl
@[simp] theorem elemPreFrobenioid_Base {A B : ElemFrobCat Φ} (f : A ⟶ B) :
    (elemPreFrobenioid Φ hD hpd).Base f = f.base := rfl

/-! ### ★`𝔽_Φ` は isotropic type

原文 (FrdI p.27):
> type]. This completes the proof of assertion (i). Assertion (iii) is immediate from
-/

/-- ★★**`𝔽_Φ` のすべての対象は isotropic**(一般)。

isometric な pre-step `φ` は
「零因子 `φ.div` が**可逆**(`toChar_eq_zero_iff`)・次数 1・底が同型」なので、
**逆射を明示的に構成できる**: `⟨inv φ.base, Φ.map (inv φ.base) U.neg, 1⟩`。

★`wP` では「零因子が **0**」だったが、一般には「**可逆**」である。
`Φ = ℕ` が sharp なので両者が一致していた——**これが模型固有の性質だった**。 -/
theorem elemFrob_isotropic (A : ElemFrobCat Φ) :
    IsIsotropic (elemPreFrobenioid Φ hD hpd) A := by
  intro Dd φ hiso hstep
  have hdeg : φ.deg = 1 := hstep.1
  haveI hbase : IsIso φ.base := hstep.2
  have hu : IsAddUnit φ.div := toChar_eq_zero_iff.mp hiso
  obtain ⟨U, hU⟩ := hu
  refine ⟨⟨inv φ.base, Φ.map (inv φ.base) U.neg, 1⟩, ?_, ?_⟩
  · refine ElemFrobCat.Hom.ext (IsIso.hom_inv_id _) ?_ ?_
    · show Φ.map φ.base (Φ.map (inv φ.base) U.neg) + ((1 : ℕ+) : ℕ) • φ.div
        = (0 : Φ.val A.base)
      rw [← Φ.map_comp, IsIso.hom_inv_id, Φ.map_id]
      simp only [PNat.one_coe, one_smul, ← hU]
      rw [add_comm]
      exact U.val_neg
    · show (1 : ℕ+) * φ.deg = 1
      simp [hdeg]
  · refine ElemFrobCat.Hom.ext (IsIso.inv_hom_id _) ?_ ?_
    · show Φ.map (inv φ.base) φ.div + (φ.deg : ℕ) • Φ.map (inv φ.base) U.neg
        = (0 : Φ.val Dd.base)
      rw [hdeg]
      simp only [PNat.one_coe, one_smul, ← map_add, ← hU]
      rw [U.val_neg, map_zero]
    · show φ.deg * 1 = 1
      simp [hdeg]

/-- ★**`𝔽_Φ` のすべての射は co-angular**。

★これは `Proposition 1.4, (i)` の直接の帰結であり、**原文自身が同じ道を通る**。

原文 (FrdI p.27):
> sition 1.4, (i) [which is applicable to all morphisms of FΦ since FΦ is of isotropic
-/
theorem elemFrob_coAngular {A B : ElemFrobCat Φ} (φ : A ⟶ B) :
    IsCoAngular (elemPreFrobenioid Φ hD hpd) φ :=
  prop_1_4_i _ φ (fun X _ => elemFrob_isotropic Φ hD hpd X)

/-- ★**`Proposition 1.4, (i)` の「In particular」を `𝔽_Φ` に当てたもの** ——
Frobenius 型 ⟺ isometric な base-isomorphism。 -/
theorem elemFrob_frobeniusType_iff {A B : ElemFrobCat Φ} (φ : A ⟶ B) :
    IsFrobeniusType (elemPreFrobenioid Φ hD hpd) φ ↔
      (IsIsometric (elemPreFrobenioid Φ hD hpd) φ ∧
        IsBaseIsomorphism (elemPreFrobenioid Φ hD hpd) φ) :=
  prop_1_4_i_frobeniusType _ φ (fun X _ => elemFrob_isotropic Φ hD hpd X)

/-! ### ★`.src` を**あえて付けない**

★`Proposition 1.5` は **まだ完成していない**——(i) の `Definition 1.3` 22条、
(ii) の `𝒪^▷(A) ≅ Φ(A)`、7つの type、(iii) がすべて未了である。

★`tools/frdi-progress.mjs` は「`.src` が項目を指していれば実装済み」と数えるので、
ここで `.src` を付けると **`Proposition 1.5` が完了として数えられてしまう**。
**部分実装を完了として数えさせない**ために、あえて付けない。

★同じ判断を `MonoidOn.charOn` でもした。**やった仕事が数に出ない**方が、
**やっていない仕事が数に出る**より良い。
-/

end ABC3.Found.FrdI
