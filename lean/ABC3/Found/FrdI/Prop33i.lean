import ABC3.Found.FrdI.Remark312
import ABC3.Found.FrdI.Prop22

/-!
# [FrdI] Proposition 3.3, (i) —— `End(𝒞^pl-bk_A → 𝒞)^bs-iso`

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.59–p.61。

原文 (FrdI p.59):
> Proposition 3.3. (Base-identity Pre-steps and Units)

## ★★(i) の主張は 2 つ

原文 (FrdI p.60):
> Proof. First, we consider assertion (i). Note that since the composite of the functor

| # | 内容 |
|---|---|
| 1 | `𝒟` が Frobenius-slim なら、`𝔽 → End(𝒞^pl-bk_A → 𝒞)^bs-iso` による `1 ∈ ℤ≥0` の像の成分は**すべて base-identity pre-step** |
| 2 | 逆に `𝒞` が Frobenius-normalized 型で `A` が Frobenius-trivial なら、`A` のどの base-identity pre-step 自己射もそう現れる |

## ★★★主張 1 の骨(原文どおり、2026-08-17 に測定)

原文 (FrdI p.60):
> obtained by considering the Frobenius degree of the induced endomorphism of A

★**2 つの準同型に分けて、それぞれ `𝔽 ↠ ℕ≥1` を経由させる**:

| 成分 | 行き先 | 経由の理由 |
|---|---|---|
| **底** | `Aut(𝒟_{A_𝒟} → 𝒟)` | ★`𝒟` が **Frobenius-slim**(定義そのもの) |
| **次数** | `ℕ≥1` | ★★**`ℕ≥1` が可換かつ簡約的**(`elemFrob_factors_of_cancel`) |

★どちらも `⟨1,1⟩` の次数が `1` なので、像は `1` になる。
★**底が恒等 = base-identity、次数が 1 = linear。合わせて base-identity pre-step**
(すなわち `𝒪^▷(−)` の元)である。

★★底の側で `𝒞^pl-bk_A → 𝒟` が `𝒟_{A_𝒟} → 𝒟` を経由することは
**`Definition 1.3, (i), (c)` の圏同値**(`plBkEquiv`)である。

## ★本ファイルの範囲

★**土台(対象の定義)と、次数の側**をここで取る。
★底の側は `plBkEquiv` を通した `Aut` の移送が要る(別途)。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

/-! ## ★自然な関手 `𝒞^pl-bk_A → 𝒞` -/

/-- ★**`𝒞^pl-bk_A → 𝒞`** —— pull-back 射の広い部分圏のスライスから `𝒞` へ戻る。

★`map_id` も `map_comp` も **`rfl`** である(広い部分圏は射を包んでいるだけだから)。 -/
def plBkOverToC (A : C) : Over (⟨A⟩ : PlBk P) ⥤ C where
  obj Z := Z.left.obj
  map f := f.left.hom
  map_id _ := rfl
  map_comp _ _ := rfl

/-- ★**底へ落とすと `plBkOverFunctor` を経由する** —— `rfl`。

★★これが原文の「the composite of the functor `𝒞^pl-bk_A → 𝒞` with the natural
projection functor `𝒞 → 𝒟` factors as …」の中身である。 -/
theorem plBkOverToC_comp_proj (A : C) :
    plBkOverToC P A ⋙ P.proj = plBkOverFunctor P A ⋙ Over.forget _ := rfl

/-! ## ★`End(𝒞^pl-bk_A → 𝒞)^bs-iso` -/

/-- ★★**原文の `End(𝒞^pl-bk_A → 𝒞)^bs-iso`** ——
成分がすべて base-isomorphism である自然変換のなす部分モノイド。

★`End F` の積は `x * y = y ≫ x` なので、`mul_mem` は
**base-isomorphism が合成で閉じる**ことに帰着する。 -/
def endBsIso (A : C) : Submonoid (End (plBkOverToC P A)) where
  carrier := {η | ∀ Z, IsBaseIsomorphism P (NatTrans.app η Z)}
  one_mem' := fun Z => by
    show IsIso (P.Base (NatTrans.app (𝟙 (plBkOverToC P A)) Z))
    rw [show NatTrans.app (𝟙 (plBkOverToC P A)) Z = 𝟙 _ from rfl, P.Base_id]
    infer_instance
  mul_mem' {x y} hx hy := fun Z => by
    show IsIso (P.Base (NatTrans.app (y ≫ x) Z))
    rw [show NatTrans.app (y ≫ x) Z = NatTrans.app y Z ≫ NatTrans.app x Z from rfl]
    exact isBaseIsomorphism_comp P (hy Z) (hx Z)

/-! ## ★次数の側 —— `End(…)^bs-iso →* ℕ≥1`

★**恒等対象 `(A ⟶ A)`(恒等射は pull-back)での成分**の Frobenius 次数を取る。
★次数は乗法的(`degFr_comp`)なので準同型になる。 -/

/-- ★恒等射が定める `𝒞^pl-bk_A` の対象。 -/
def plBkIdObj (A : C) : Over (⟨A⟩ : PlBk P) :=
  Over.mk (⟨𝟙 A, isPullBack_of_isIso P (𝟙 A)⟩ :
    (⟨A⟩ : PlBk P) ⟶ (⟨A⟩ : PlBk P))

/-- ★★**次数を取る準同型** `End(𝒞^pl-bk_A → 𝒞) →* ℕ≥1`。

★`End` の積が `x * y = y ≫ x` で、`degFr` が `degFr (ψ ≫ φ) = degFr φ * degFr ψ`
なので、**両方の反転が打ち消して**準同型になる。 -/
def endDegHom (A : C) : End (plBkOverToC P A) →* ℕ+ where
  toFun η := P.degFr (NatTrans.app η (plBkIdObj P A))
  map_one' := by
    show P.degFr (NatTrans.app (𝟙 (plBkOverToC P A)) (plBkIdObj P A)) = 1
    rw [show NatTrans.app (𝟙 (plBkOverToC P A)) (plBkIdObj P A) = 𝟙 _ from rfl]
    exact degFr_of_isIso P (𝟙 _)
  map_mul' x y := by
    show P.degFr (NatTrans.app (y ≫ x) (plBkIdObj P A)) = _
    rw [show NatTrans.app (y ≫ x) (plBkIdObj P A)
      = NatTrans.app y (plBkIdObj P A) ≫ NatTrans.app x (plBkIdObj P A) from rfl,
      P.degFr_comp]

/-- ★★★**主張 1 の「次数の側」** —— `𝔽 → End(𝒞^pl-bk_A → 𝒞)` による
`⟨1,1⟩` の像は **linear**(次数 1)である。

原文 (FrdI p.60):
> obtained by considering the Frobenius degree of the induced endomorphism of A

★★**`ℕ≥1` が可換かつ簡約的**なので `elemFrob_factors_of_cancel` が効き、
合成 `𝔽 → ℕ≥1` は `deg` を経由する。★`⟨1,1⟩` の次数は `1` なので像も `1`。 -/
theorem endDeg_gen_eq_one (A : C) (f : ElemFrob ℕ →* End (plBkOverToC P A)) :
    P.degFr (NatTrans.app (f ⟨1, 1⟩) (plBkIdObj P A)) = 1 := by
  obtain ⟨g, hg⟩ := elemFrob_factors_of_cancel ((endDegHom P A).comp f)
  have h := congrArg (fun m : ElemFrob ℕ →* ℕ+ => m ⟨1, 1⟩) hg
  show (endDegHom P A) (f ⟨1, 1⟩) = 1
  show ((endDegHom P A).comp f) ⟨1, 1⟩ = 1
  rw [h]
  show g ((⟨1, 1⟩ : ElemFrob ℕ).deg) = 1
  show g 1 = 1
  exact g.map_one

end ABC3.Found.FrdI
