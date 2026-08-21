/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.TypeTransport
import ABC3.Found.FrdI.Cor54Rigid
import ABC3.Found.FrdI.Cor57Model
import ABC3.Found.FrdI.Prop32Perfect
import ABC3.Found.FrdI.MonoidTransport
import ABC3.Found.FrdI.Prop55UntrCat

/-!
# [FrdI] Proposition 5.5, (iii) —— standard 型が `𝒞^un-tr` へ移る

原文 (FrdI p.105):
> the pull-back morphisms of Cun-tr, Crlf are precisely the linear isometries [cf. Proposition 1.4,

★★原文 (iii) の最後の文は
「`𝒞` が not group-like で standard(resp. rationally standard)型なら
`𝒞^un-tr`・`𝒞^rlf` もそう」である。

`IsOfStandardType` は 6 条:

| 条 | `𝒞^un-tr` での根拠 |
|---|---|
| quasi-isotropic | `unTr_isotropic`(すべての対象が isotropic) |
| Frobenius-isotropic | 恒等射が Frobenius 型 ＋ isotropic |
| **Frobenius-normalized** | ★本ファイル(`unTr_frobNormalizedType`) |
| `𝒟` が FSMFF | `𝒞` と同じ `𝒟` |
| `Φ` が non-dilating | `𝒞` と同じ `Φ` |
| group-like ⟹ compact 対象 | ★未(`𝒞` が not group-like なら前件が偽になることを要する) |

## ★★Frobenius-normalized の筋

★`𝒞^un-tr` は **model 型**(`unTr_isOfModelType`)なので、
`Theorem 5.2, (iv)` の圏同値 `unTr_modelFrobenioid` で model Frobenioid と結べる。
★model Frobenioid では Frobenius-normalized は**無条件**である
(`ModelData.model_frobNormalizedType` —— 射が 4 成分の明示的な組だから)。

★★★そこで**逆向きに引き戻す**。要る道具は `TypeTransport.lean` の

  `isFrobeniusNormalized_map_of_toElem_iso`

で、仮定は「充満忠実 ＋ `𝔽_Φ` への関手と 1-可換」だけである。
1-可換性は在庫の `modelTypeEquiv_comp_toElem`(`Cor54Rigid.lean`)がそのまま与える。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2 u4 v4

/-! ## ★1. `RatFnData` からの `Hyp` と、1-可換性の逆向き -/

section Generic

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}

/-- ★★`RatFnData` からは `ModelData.Hyp` が**無料で出る**。

★`B` が group-like であることは `RatFnData` の `bneg` がそのまま与え、
残る 3 条(`Φ` が divisorial・`𝒟` が totally epimorphic・`𝒟` が connected)は
`PreFrobenioid` のフィールドそのものである。 -/
theorem RatFnData.hyp {G : Frobenioid P} (R : RatFnData P G) : ModelData.Hyp R.model where
  divisorial := P.divisorial
  bmonGroupLike A := (isGroupLike_iff _).mpr fun a =>
    ⟨⟨a, R.bneg A a, (add_comm a (R.bneg A a)).trans (R.bneg_add A a), R.bneg_add A a⟩, rfl⟩
  totEpiD := P.totEpiD
  connectedD := P.connectedD

end Generic

section InverseIso

variable {D₀ : Type u} [Category.{v} D₀] {Φ₀ : MonoidOn.{v, u, w} D₀}
  {C₁ : Type u2} [Category.{v2} C₁] {C₂ : Type u4} [Category.{v4} C₂]
  {P₁ : PreFrobenioid C₁ Φ₀} {P₂ : PreFrobenioid C₂ Φ₀}

/-- ★★**`𝔽_Φ` への関手との 1-可換性は圏同値の逆向きへ移る**。

★`η` を `e.inverse` で左から whisker し、counit で `𝟭` に潰すだけ。 -/
noncomputable def inverse_comp_toElem_iso (e : C₁ ≌ C₂)
    (η : e.functor ⋙ P₂.toElem ≅ P₁.toElem) :
    e.inverse ⋙ P₁.toElem ≅ P₂.toElem :=
  Functor.isoWhiskerLeft e.inverse η.symm ≪≫ (Functor.associator _ _ _).symm
    ≪≫ Functor.isoWhiskerRight e.counitIso P₂.toElem ≪≫ Functor.leftUnitor _

end InverseIso

/-! ## ★2. `𝒞^un-tr` は Frobenius-normalized 型 -/

section UnTr

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} [IsConnected D]

set_option maxHeartbeats 800000 in
/-- ★★★★**[FrdI] Proposition 5.5, (iii) の一条** —— `𝒞^un-tr` は **Frobenius-normalized 型**。

★★筋は 3 段:
1. `𝒞^un-tr` は model 型(`unTr_isOfModelType`)なので model Frobenioid と圏同値。
2. model Frobenioid では Frobenius-normalized は**無条件**
   (`ModelData.model_frobNormalizedType`)。
3. 1-可換性(`modelTypeEquiv_comp_toElem`)を**逆向き**にして引き戻す
   (`isFrobeniusNormalized_map_of_toElem_iso`)。 -/
theorem unTr_frobNormalizedType (Fc : FrobenioidCore P) (G : Frobenioid P)
    (hint : ∀ A : D, IsIntegralMonoid (Φ.val A)) (hfsmD : IsOfFSMType D) :
    IsOfFrobeniusNormalizedType (unTrPre P Fc) := by
  intro A
  have h : ModelData.Hyp (unTr_ratFnData Fc G hint hfsmD).model :=
    RatFnData.hyp (unTr_ratFnData Fc G hint hfsmD)
  have η : (unTr_modelFrobenioid Fc G hint hfsmD).functor
      ⋙ (ModelData.modelPre h).toElem ≅ (unTrPre P Fc).toElem :=
    modelTypeEquiv_comp_toElem (unTr_ratFnData Fc G hint hfsmD)
      (unTr_isotropic P Fc) (unTr_isOfModelType Fc G)
  have η' := inverse_comp_toElem_iso (unTr_modelFrobenioid Fc G hint hfsmD) η
  have hback : IsFrobeniusNormalized (unTrPre P Fc)
      ((unTr_modelFrobenioid Fc G hint hfsmD).inverse.obj
        ((unTr_modelFrobenioid Fc G hint hfsmD).functor.obj A)) :=
    isFrobeniusNormalized_map_of_toElem_iso _ η' _
      (ModelData.model_frobNormalizedType h _)
  exact isFrobeniusNormalized_of_iso
    ((unTr_modelFrobenioid Fc G hint hfsmD).unitIso.app A).symm hback

/-- ★★`𝒞^un-tr` は **Frobenius-isotropic 型**(すべての対象が isotropic だから恒等射でよい)。 -/
theorem unTr_frobIsotropicType (Fc : FrobenioidCore P) :
    IsOfFrobeniusIsotropicType (unTrPre P Fc) :=
  fun A => ⟨A, 𝟙 A, isFrobeniusType_of_isIso (unTrPre P Fc) (𝟙 A), unTr_isotropic P Fc A⟩

/-! ## ★3. standard 型が `𝒞^un-tr` へ移る

★★`IsOfStandardType` の 6 条のうち 5 条は上と在庫で埋まり、残る
`groupLikeCompact` は**前件が偽**になる ——
原文の「suppose that `𝒞`, hence also `𝒞^un-tr`, `𝒞^rlf`, are not of group-like type」が
それである。

★`𝒞^un-tr` の対象は `𝒞^istr` の対象そのもの(`UnTr P := Istr P`)で、
`𝔽_Φ` への関手も `A ↦ P.toElem.obj A.obj` なので、
**group-like 性は `𝒞` のそれと同じ条件**である(`IsGroupLikeObj` は `Φ` と底だけで決まる)。
★したがって `𝒞` が group-like 型なら `𝒞^un-tr` もそう(`unTr_isOfGroupLikeType`)。
★★逆(原文の「hence also」)は `𝒞` の各対象を isotropic 包へ送る段が要るので、
ここでは **`¬ IsOfGroupLikeType (unTrPre P Fc)` を仮引数で受ける**。 -/

omit [IsConnected D] in
/-- ★**`𝒞` が group-like 型なら `𝒞^un-tr` も** —— 対象も底も同じだから。 -/
theorem unTr_isOfGroupLikeType (Fc : FrobenioidCore P) (h : IsOfGroupLikeType P) :
    IsOfGroupLikeType (unTrPre P Fc) :=
  fun A => h (show Istr P from A).obj

set_option maxHeartbeats 800000 in
/-- ★★★★★**[FrdI] Proposition 5.5, (iii)** —— **standard 型は `𝒞^un-tr` へ移る**。

原文 (FrdI p.105):
> if C is of standard (respectively, rationally standard) type, then so are Cun-tr, Crlf.

★★6 条の内訳:
| 条 | 根拠 |
|---|---|
| quasi-isotropic | `unTr_isotropic`(isotropic 型 ⟹ quasi-isotropic) |
| Frobenius-isotropic | `unTr_frobIsotropicType` |
| group-like ⟹ compact 対象 | ★**前件が偽**(`hngl`) |
| Frobenius-normalized | `unTr_frobNormalizedType` |
| `𝒟` が FSMFF | `𝒞` と同じ `𝒟` |
| `Φ` が non-dilating | `𝒞` と同じ `Φ` | -/
theorem unTr_standardType (Fc : FrobenioidCore P) (G : Frobenioid P)
    (hint : ∀ A : D, IsIntegralMonoid (Φ.val A)) (hfsmD : IsOfFSMType D)
    (F' : FrobenioidCore (unTrPre P Fc))
    (hngl : ¬ IsOfGroupLikeType (unTrPre P Fc))
    (hstd : IsOfStandardType D C P Fc) :
    IsOfStandardType D (UnTr P) (unTrPre P Fc) F' where
  quasiIsotropic :=
    isOfQuasiIsotropicType_of_isOfIsotropicType (unTrPre P Fc) F' (unTr_isotropic P Fc)
  frobIsotropic := unTr_frobIsotropicType Fc
  groupLikeCompact hg := absurd hg hngl
  frobNormalized := unTr_frobNormalizedType Fc G hint hfsmD
  baseFSMFF := hstd.baseFSMFF
  phiNonDilating := hstd.phiNonDilating

end UnTr

/-! ## ★4. standard 型が `𝒞^rlf` へ移る

★★`𝒞^rlf` は **model Frobenioid そのもの**(`ScModelObj S G … = ModelData.Obj (scModel …)`)
なので、`ModelData.model_isOfStandardType` が**直接当たる** —— 移送は要らない。
★残る入力は 3 つ:
* `hfsmff` —— `𝒟` が FSMFF(`𝒞` と同じ `𝒟` なので `𝒞` の standard 性から)
* `hnd` —— **係数拡大した `Φ` が non-dilating**(★実化が non-dilating を保つこと。未)
* `hngl` —— `𝒞^rlf` が not group-like(原文の「hence also」。未)
-/

section Rlf

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {S : Type} [CommSemiring S]

set_option maxHeartbeats 800000 in
/-- ★★★★★**[FrdI] Proposition 5.5, (iii)** —— **standard 型は `𝒞^rlf` へ移る**。

★`𝒞^rlf` は model Frobenioid そのものなので、`model_isOfStandardType` の直接適用である
(quasi-isotropic・Frobenius-isotropic・Frobenius-normalized の 3 条は
model Frobenioid では**無条件**)。 -/
theorem scModel_standardType (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    (hcharInj : ∀ {A B : D} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S) (Φ.map α)))
    (hint : ∀ A : D, IsIntegralMonoid (ScT S (Φ.val A)))
    (hfsmD : IsOfFSMType D)
    (hdiv : (scModel S G hiso hfn hcharInj hint hfsmD).phi.IsDivisorialOn)
    (htot : IsTotallyEpimorphic D) (hconn : IsConnected D)
    (F' : FrobenioidCore (ModelData.modelPre
      (scModel_hyp G hiso hfn hcharInj hint hfsmD hdiv htot hconn)))
    (hfsmff : IsOfFSMFFType D)
    (hnd : (scModel S G hiso hfn hcharInj hint hfsmD).phi.IsNonDilatingOn)
    (hngl : ¬ IsOfGroupLikeType (ModelData.modelPre
      (scModel_hyp G hiso hfn hcharInj hint hfsmD hdiv htot hconn))) :
    IsOfStandardType D (ScModelObj S G hiso hfn hcharInj hint hfsmD)
      (ModelData.modelPre (scModel_hyp G hiso hfn hcharInj hint hfsmD hdiv htot hconn)) F' :=
  ModelData.model_isOfStandardType _ F' hfsmff hnd (fun hg => absurd hg hngl)

end Rlf

/-! ## ★5. standard 型が `𝒞^pf` へ移る

原文 (FrdI p.105):
> assertion (iii) for Cpf follows immediately from the deﬁnitions [cf. also Proposition

★★原文が名指しするのは「`𝒞^pf` が isotropic 型であること」と
「assertion (i) により `𝒪^▷` の像が `𝒪^▷` の perfection であること」の 2 つである。
★前者は在庫(`pfRoot_isOfIsotropicType`)にあり、`IsOfStandardType` の 6 条のうち
**4 条が無料で埋まる**。

★★残る 3 条は葉として名前をつけておく:

1. ★**`Φ^pf` が non-dilating** —— ★★**閉じた**(`MonoidOn.pfOn_isNonDilatingOn`、下の節)。
   ★★測って分かったこと: `MPrec a b` は `a ≤ n·b`(**正の倍数まで許す**)なので、
   `M^pf` の元 `a` と `a/2` は `MPrec` について**同値**である。
   したがって「perfect だから primary 元が消える」ということは**ない**
   (一度そう読み違えたので記す)。★筋は
   `(M^pf)^char ≅ (M^char)^pf` と「`M^char → (M^pf)^char` が `MPrec` を保ち反射する」で、
   `PfGp.lean` の `(M^pf)^gp ≅ (M^gp)^pf` と同じ形の仕事になる。
2. **`𝒞^pf` が Frobenius-normalized 型**(`hfn`)。
   ★`𝒞` の Frobenius-normalized 性を代表元に降ろし、共通の添字で当てて押し戻す段。
3. **group-like のときの Frobenius-compact 対象**(`hgl`)。
   ★`𝒞` が not group-like なら前件が偽になる(`𝒞^un-tr` の側と同じ形)。 -/

/-! ### ★5-a. perfection は non-dilating 性を保つ

★★これで上の 3 条のうち **1 番目が閉じる**。

★★★筋は 3 段:
1. `M^pf` は sharp(`Pf.isSharp_pf`)なので、non-dilating は
   **`M^pf` の上の含意**に降りる(`isNonDilating_of_sharp`)。
2. `Pf.of` は準素元を準素元へ送り(`Pf.of_primary`)、`≼` を**保ち反射する**
   (`Pf.of_mprec` / `Pf.mprec_of_of`)。★反射は在庫の `Pf.mle_num_of_mle`
   (「`Pf` の `≼` は分子を `k` 倍すれば `M` の `≼` に落ちる」)＋
   「`x ≼ k·x`」で出る。
3. したがって `M^pf` 側の仮定が `M` 側の仮定に落ち、`M` の non-dilating が
   `α = id` を与え、`Pf.map_id` で戻る。 -/

section PfMonoid

variable {M : Type w} [AddCommMonoid M]

theorem mLe_trans' {a b c : M} (h₁ : MLe a b) (h₂ : MLe b c) : MLe a c := by
  obtain ⟨x, hx⟩ := h₁
  obtain ⟨y, hy⟩ := h₂
  exact ⟨x + y, by rw [← add_assoc, hx, hy]⟩

theorem mLe_nsmul_self (x : M) (k : ℕ) (hk : 0 < k) : MLe x (k • x) := by
  refine ⟨(k - 1) • x, ?_⟩
  rw [← succ_nsmul' x (k - 1)]
  congr 1
  omega

theorem mLe_nsmul {x y : M} (k : ℕ) (h : MLe x y) : MLe (k • x) (k • y) := by
  obtain ⟨c, hc⟩ := h
  exact ⟨k • c, by rw [← smul_add, hc]⟩

theorem mPrec_of_nsmul {x y : M} {k : ℕ} (hk : 0 < k) (h : MPrec (k • x) y) : MPrec x y := by
  obtain ⟨n, hn, hle⟩ := h
  exact ⟨n, hn, mLe_trans' (mLe_nsmul_self x k hk) hle⟩

theorem mPrec_nsmul_right {x y : M} {k : ℕ} (hk : 0 < k) (h : MPrec x (k • y)) :
    MPrec x y := by
  obtain ⟨n, hn, hle⟩ := h
  exact ⟨n * k, Nat.mul_pos hn hk, by rwa [mul_smul]⟩

theorem Pf.of_eq_nsmul_mk (x : M) (k : ℕ+) :
    Pf.of x = ((k : ℕ+) : ℕ) • Pf.mk x k := by
  rw [Pf.nsmul_mk]
  refine (Pf.sound 1 ?_).symm
  push_cast
  simp [smul_smul]

/-- ★`Pf.of` は `≼` を保つ。 -/
theorem Pf.of_mprec {x y : M} (h : MPrec x y) : MPrec (Pf.of x) (Pf.of y) := by
  obtain ⟨n, hn, c, hc⟩ := h
  refine ⟨n, hn, Pf.of c, ?_⟩
  rw [← map_add, hc, map_nsmul]

/-- ★★`Pf.of` は `≼` を**反射する** —— 在庫の `Pf.mle_num_of_mle` がそのまま効く。 -/
theorem Pf.mprec_of_of {x y : M} (h : MPrec (Pf.of x) (Pf.of y)) : MPrec x y := by
  obtain ⟨n, hn, hle⟩ := h
  have hle' : MLe (Pf.mk x 1) (Pf.mk (n • y) 1) := by simpa using hle
  obtain ⟨k, c, hc⟩ := Pf.mle_num_of_mle hle'
  refine mPrec_of_nsmul (k := ((k : ℕ+) : ℕ)) k.pos
    ⟨((k : ℕ+) : ℕ) * n, Nat.mul_pos k.pos hn, c, ?_⟩
  rw [hc, mul_smul]

/-- ★★★**`Pf.of` は準素元を準素元へ送る**(`M` が divisorial なら)。

★★分母は `≼` にとって**見えない** —— `Pf.of c₀ = k • (c₀/k)` なので、
`c₀/k` と `c₀` は `≼` について同値である。 -/
theorem Pf.of_primary (hdiv : IsDivisorial M) {b : M} (h : IsPrimaryElt b) :
    IsPrimaryElt (Pf.of b) := by
  refine ⟨?_, ?_⟩
  · intro h0
    refine h.1 ((Pf.mk_eq_mk_same_iff hdiv (r := 1)).mp ?_)
    show Pf.mk b 1 = Pf.mk 0 1
    exact h0
  · intro c hc hcb
    induction c using Pf.inductionOn with | _ c₀ k =>
    have hc₀ : c₀ ≠ 0 := by
      intro h0
      exact hc (by rw [h0]; exact Pf.mk_zero k)
    have h1 : MPrec (Pf.of c₀) (Pf.of b) := by
      obtain ⟨n, hn, hle⟩ := hcb
      refine ⟨((k : ℕ+) : ℕ) * n, Nat.mul_pos k.pos hn, ?_⟩
      rw [Pf.of_eq_nsmul_mk c₀ k, mul_smul]
      exact mLe_nsmul _ hle
    have h4 : MPrec (Pf.of b) (Pf.of c₀) :=
      Pf.of_mprec (h.2 c₀ hc₀ (Pf.mprec_of_of h1))
    refine mPrec_nsmul_right (k := ((k : ℕ+) : ℕ)) k.pos ?_
    rwa [← Pf.of_eq_nsmul_mk c₀ k]

/-- ★★★★★**perfection は non-dilating 性を保つ**。 -/
theorem Pf.isNonDilating_map (hdiv : IsDivisorial M) (α : M →+ M) (hα : IsNonDilating α) :
    IsNonDilating (Pf.map α) := by
  refine isNonDilating_of_sharp (Pf.isSharp_pf hdiv.2) (Pf.map α) (fun h => ?_)
  have hM : ∀ b : M, IsPrimaryElt b → MPrec (α b) b := by
    intro b hb
    have hkey := h (Pf.of b) (Pf.of_primary hdiv hb)
    refine Pf.mprec_of_of ?_
    have he : Pf.map α (Pf.of b) = Pf.of (α b) := Pf.map_mk α b 1
    rwa [he] at hkey
  rw [addMonoidHom_eq_id_of_primary_mprec hdiv.2 α hα hM, Pf.map_id]

end PfMonoid

section PfOn

variable {D : Type u} [Category.{v} D]

/-- ★★★★★**`Φ^pf` は non-dilating**(`Φ` が non-dilating で各値が divisorial なら)。 -/
theorem MonoidOn.pfOn_isNonDilatingOn (Φ : MonoidOn.{v, u, w} D)
    (hsh : ∀ A : D, IsSharp (Φ.val A)) (hdiv : ∀ A : D, IsDivisorial (Φ.val A))
    (h : Φ.IsNonDilatingOn) : (Φ.pfOn hsh).IsNonDilatingOn := by
  intro A e
  show IsNonDilating (Pf.map (Φ.map e))
  exact Pf.isNonDilating_map (hdiv A) _ (h A e)

end PfOn

section Pf

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {F : FrobenioidCore P}

set_option maxHeartbeats 800000 in
/-- ★★★★★**[FrdI] Proposition 5.5, (iii)** —— **standard 型は `𝒞^pf` へ移る**。

★6 条のうち 4 条(quasi-isotropic・Frobenius-isotropic・`𝒟` が FSMFF)は
`pfRoot_isOfIsotropicType` と「`𝒟` が同じ」で無料。
★残る 2 条(Frobenius-normalized・group-like のときの compact 対象)は
仮引数で受けている(上の節に筋を記した)。 -/
theorem pfRoot_standardType (hfi : IsOfFrobeniusIsotropicType P)
    (F₂ : FrobenioidCore (pfRootPre P F))
    (hfn : IsOfFrobeniusNormalizedType (pfRootPre P F))
    (hgl : IsOfGroupLikeType (pfRootPre P F) →
      ∃ A : Istr (pfRootPre P F), IsFrobeniusCompact (istrPre (pfRootPre P F) F₂) A)
    (hstd : IsOfStandardType D C P F) :
    IsOfStandardType D (PfRootObj P F) (pfRootPre P F) F₂ where
  quasiIsotropic :=
    isOfQuasiIsotropicType_of_isOfIsotropicType (pfRootPre P F) F₂
      (pfRoot_isOfIsotropicType (F := F) hfi)
  frobIsotropic := fun A => ⟨A, 𝟙 A, isFrobeniusType_of_isIso (pfRootPre P F) (𝟙 A),
    pfRoot_isOfIsotropicType (F := F) hfi A⟩
  groupLikeCompact := hgl
  frobNormalized := hfn
  baseFSMFF := hstd.baseFSMFF
  phiNonDilating :=
    MonoidOn.pfOn_isNonDilatingOn Φ (phiSharp P) P.divisorial hstd.phiNonDilating

end Pf

/-! ### ★出典の紐付け -/

/-- ★★★★locator —— `Proposition 5.5, (iii)` の「`𝒞^un-tr` は Frobenius-normalized 型」の条。 -/
def unTr_frobNormalizedType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (iii) — 𝒞^un-tr は Frobenius-normalized 型",
    sectionId := "frdi-prop-5-5" }

/-- ★★★★★locator —— `Proposition 5.5, (iii)` の「standard 型が `𝒞^un-tr` へ移る」の条
(★**条つき**: `𝒞^un-tr` が not group-like であることを仮引数で受けている)。 -/
def unTr_standardType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (iii) — standard 型は 𝒞^un-tr へ移る",
    sectionId := "frdi-prop-5-5" }

/-- ★★★★★locator —— `Proposition 5.5, (iii)` の「standard 型が `𝒞^pf` へ移る」の条
(★**条つき**: `Φ^pf` の non-dilating・Frobenius-normalized・group-like のときの
compact 対象の 3 条を仮引数で受けている)。 -/
def pfRoot_standardType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (iii) — standard 型は 𝒞^pf へ移る",
    sectionId := "frdi-prop-5-5" }

end ABC3.Found.FrdI
