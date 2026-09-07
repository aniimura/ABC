import ABC3.Found.PGC.ArtinMap

/-!
# 経路 Λ9 —— `tors_{p^n}(Gal(K^ab/K)) ≅ μ_{p^n}(K)`

典拠: S. Mochizuki, *A Version of the Grothendieck Conjecture for p-adic Local Fields*
(1997) Section 1、Proposition 1.1 の直後の地の文（物理 p.3）。構造化済み原文は
`ResearchPaper/1_Structured/A Version of the Grothendieck Conjecture for p-adic Local Fields/section-1.html`
の `#prop-1-1`。

原文 (pGC p.3):
> Now recall from local class field theory (see, e.g., [3]) that we have a natural isomorphism ΓKab ≅ (K×)∧ where the superscripted "ab" denotes the abelianization, and "(K×)∧" denotes the profinite completion of K× = K − {0}. Let k be the residue field of OK (the ring of integers of K). Thus, k is the field of q = pf elements. Now it is well-known that (K×)∧ fits into an exact sequence of topological groups: 0 → UK → (K×)∧ → Z → 0 where UK d=ef OK×.

## 何を作ったか

pGC Proposition 1.1（円分指標の群論的復元）が名指しで推奨する道の最後の環である。

| 段 | 宣言 | 中身 |
|---|---|---|
| 1 | `powTorsionProdEquiv` | ★抽象核（純群論）: `B` が `m` 捩れ自由なら `tors_m(A × B) ≅ tors_m A` |
| 2 | `ContinuousMulEquiv.prodCongr` | ★抽象核（位相）: 位相群同型の直積 |
| 3 | `galEquivUnitsZHat` / `galContinuousEquivUnitsZHat` | `Gal(K^ab/K) ≅ 𝒪_K^× × Ẑ`（位相群として） |
| 4 | `torsionEquivRootsOfUnity` | `tors_m G ≅ μ_m(K)`（`G ≅ 𝒪_K^× × Ẑ` を持つ任意の可換群 `G`） |
| 5 | `exists_abelianGalTorsion_equiv_rootsOfUnity` | ★★到達点（仮説なし） |
| 6 | `torsionEquivRootsOfUnity_cyclotomic` | 作用が円分指標倍であること |
| 7 | `natCard_powTorsion_abelianGal_dvd` 他 | 退化の自己検査 |

## 上流の在庫（本ファイルは新規の数学をほとんど足していない）

* **Λ8**（`ArtinMap.lean`）: `abelianGalEquivProd : Gal(K^ab/K) ≅ 𝒪_K^× × Gal(K^ur/K)`
  とその連続性 `continuous_abelianGalEquivProd` / `continuous_abelianGalEquivProd_symm`。
  ★標的が `K^ab` になるのは LKW（Yoshida08 Theorem 6.15、`LocalClassFieldTheory.lean`）による。
* `ArithFrobeniusTopGen.lean::zhatMulEquivUnramGalArith : Ẑ ≃ₜ* Gal(K^ur/K)`。
* `ProfiniteUnitsTorsion.lean::torsionUnitsZHatEquiv : tors_m(𝒪_K^× × Ẑ) ≅ μ_m(K)`
  と `zhat_eq_one_of_pow_eq_one`（`Ẑ` は捩れ自由）。

★★**「`Gal(K^ab/K) ≃ₜ* 𝒪^× × Ẑ` は 3 行で出る」という Λ8 の見立ての実測**:
群同型（`galEquivUnitsZHat`）は**本当に 3 行**だった。位相群同型（`≃ₜ*`）は
`ContinuousMulEquiv.prodCongr` が **mathlib に無かった**ため、その補題（8 行）と
仮説を落とす段（`exists_abelianGalContinuousEquivUnitsZHat`、13 行）が余分に要り、
合計で 3 行では済まなかった。★見立ては群同型については当たり、位相版については外れである。

## ★★退化の自己検査

1. **`μ_{p^n}(K)` は `K` に依存する。★「常に位数 `p^n`」ではない。**
   たとえば `K = ℚ_p`（`p` 奇）なら `μ_p(ℚ_p) = 1` である
   （`ζ ≠ 1` かつ `ζ^p = 1` なら `ℚ_p(ζ)/ℚ_p` は次数 `p−1 > 1` の分岐拡大になり
   `ζ ∉ ℚ_p`）。したがって本ファイルの到達点は
   「`tors_{p^n}(Gal(K^ab/K))` が `p^n` 個ある」ではなく
   「`μ_{p^n}(K)` と**同型**である」であり、右辺は自明群でありうる。
   ★型に出ているとおり、結論は `rootsOfUnity (p ^ n) K.carrier` という
   **`K` に依存する**群であって `ZMod (p^n)` ではない。
   `natCard_powTorsion_abelianGal_dvd` は「位数は `p^n` を**割る**」しか言わない。
2. **`Ẑ` の捩れ自由性は仮定していない。** `zhat_eq_one_of_pow_eq_one`
   （`ProfiniteUnitsTorsion.lean` で実際に証明されている）を消費している。
   `Ẑ` を捩れのある副有限群に取り替えると本ファイルの主張は偽になる。
3. **`𝒪_K^×` の中の捩れと `K^×` の中の捩れが一致すること**を落としていない。
   `torsionUnitsZHatEquiv`（在庫）の中身がまさにそれで、`m` 乗根はノルム `1` なので
   自身も逆元も `𝒪_K` に入る（`mem_carrierIntegers_of_pow_eq_one`）。
   本ファイルはその同型をそのまま運んでいるので、この段を飛ばしていない。
4. **`m ≠ 0` は落とせない。** `m = 0` なら両辺とも群全体になり主張は偽。
5. **捩れ部分群そのものは `π` に依存しない。** `powTorsion Gal(K^ab/K) (p^n)` は
   `{x | x^{p^n} = 1}`（`mem_powTorsion_abelianGal`）であり、Lubin-Tate の
   素元 `π` や級数 `f` を含まない。★依存しているのは**同型 `e` の方**だけで、
   それは `∃`（あるいは `Nonempty`）の内側に閉じ込めてある。

## ★符号（Λ8 の測定を壊していない）

本ファイルは `abelianGalEquivProd` を**そのまま**通しており、
第 2 成分の Frobenius の正規化（原典 §2.2 の幾何 Frobenius `Frob_K = ϕ^{-1}`）には
一切触れていない。捩れは第 2 成分では `1` なので（`Ẑ` が捩れ自由）、
到達点は Frobenius の向きに依存しない。★`zhatMulEquivUnramGalArith` は
**算術** Frobenius で正規化された同型だが、上の理由で符号は結論に効かない。

## 逸脱（記録）

* **`Ẑ` の与え方**は `UnramifiedZhat.lean` / `ProfiniteUnitsTorsion.lean` に従う
  （`ProfiniteCompletion` の極限記述）。古典的な `Ẑ = lim ℤ/n` との一致は主張しない。
* **`μ_{p^n}(K)` の与え方**は `rootsOfUnity (p^n) K.carrier`（`Subgroup (K^×)`）。
* **`Gal(K^ab/K)` の `CommGroup` 構造**は mathlib の scoped instance
  `IsMulCommutative → CommGroup`（`Algebra/Group/Defs.lean`）で得ている
  （`AbelianClosure.lean::instIsAbelianGaloisAbelianClosure` が
  `IsAbelianGalois K.carrier (abelianClosure K)` を与えている）。
  そのため本ファイルは `open scoped IsMulCommutative` を要する。
* **原典の `(K^×)^∧`（`K^×` の副有限完備化）を作っていない。** 原典は
  `Γ_K^ab ≅ (K^×)^∧` と書くが、本木は Λ8 が `Gal(K^ab/K) ≅ 𝒪_K^× × Gal(K^ur/K)`
  を**直積の形で**持っているので、そちらを使う。
  ★どちらでも `p^n` 捩れは `μ_{p^n}(K)` になる（原典の完全列
  `0 → UK → (K^×)^∧ → Z → 0` が分裂した形が本木の直積である）。
  ★**原典の完全列が分裂することは主張していない**——本木の直積分解は
  Lubin-Tate から**独立に**得られたものである。
-/

namespace ContinuousMulEquiv

variable {A B C D : Type*} [MulOneClass A] [MulOneClass B] [MulOneClass C] [MulOneClass D]
  [TopologicalSpace A] [TopologicalSpace B] [TopologicalSpace C] [TopologicalSpace D]

/-- ★**抽象核（位相）** —— 位相群（乗法モノイド）の同型の直積。

★2026-09-07 の在庫調査で mathlib に無かった（`.cache/mathlib-index.txt` の
`ContinuousMulEquiv` の項に `prodCongr` は無く、`#158`（同名で書いて
`already been declared` を狙う）も空振りした）。分岐も付値も体も出てこない。 -/
def prodCongr (e : A ≃ₜ* C) (f : B ≃ₜ* D) : (A × B) ≃ₜ* (C × D) where
  toMulEquiv := e.toMulEquiv.prodCongr f.toMulEquiv
  continuous_toFun := by
    exact (e.continuous_toFun.comp continuous_fst).prodMk
      (f.continuous_toFun.comp continuous_snd)
  continuous_invFun := by
    exact (e.symm.continuous_toFun.comp continuous_fst).prodMk
      (f.symm.continuous_toFun.comp continuous_snd)

@[simp] theorem prodCongr_apply (e : A ≃ₜ* C) (f : B ≃ₜ* D) (x : A × B) :
    prodCongr e f x = (e x.1, f x.2) := rfl

@[simp] theorem prodCongr_symm_apply (e : A ≃ₜ* C) (f : B ≃ₜ* D) (y : C × D) :
    (prodCongr e f).symm y = (e.symm y.1, f.symm y.2) := rfl

@[simp] theorem prodCongr_toMulEquiv (e : A ≃ₜ* C) (f : B ≃ₜ* D) :
    (prodCongr e f).toMulEquiv = e.toMulEquiv.prodCongr f.toMulEquiv := rfl

end ContinuousMulEquiv

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued Classical IsMulCommutative

/-! ## §1 抽象核 A —— 純群論

★分岐・付値・Galois の語彙は 1 つも出てこない。 -/

section AbstractCore

variable {A B G : Type*} [CommGroup A] [CommGroup B] [CommGroup G]

/-- `m` 捩れ自由な群の `m` 捩れ部分群は自明。 -/
theorem powTorsion_eq_bot_of_forall_eq_one {m : ℕ} (hB : ∀ b : B, b ^ m = 1 → b = 1) :
    powTorsion B m = ⊥ :=
  eq_bot_iff.mpr fun b hb => Subgroup.mem_bot.mpr (hB b hb)

/-- 第 2 因子が `m` 捩れ自由なら、直積の `m` 捩れ元の第 2 成分は `1`。 -/
theorem snd_eq_one_of_mem_powTorsion {m : ℕ} (hB : ∀ b : B, b ^ m = 1 → b = 1)
    {x : A × B} (hx : x ∈ powTorsion (A × B) m) : x.2 = 1 :=
  hB _ (powTorsion_snd hx)

/-- ★★**抽象核** —— `B` が `m` 捩れ自由なら `tors_m (A × B) ≅ tors_m A`。

退化の自己検査: `B` の `m` 捩れ自由性は落とせない。`B = ZMod m`（乗法的に書けば
位数 `m` の巡回群）なら左辺は右辺より `m` 倍大きい。 -/
def powTorsionProdEquiv (m : ℕ) (hB : ∀ b : B, b ^ m = 1 → b = 1) :
    ↥(powTorsion (A × B) m) ≃* ↥(powTorsion A m) where
  toFun x := ⟨(x : A × B).1, powTorsion_fst x.2⟩
  invFun a := ⟨((a : A), 1), mem_powTorsion_prod a.2 (one_pow m)⟩
  left_inv x := Subtype.ext (Prod.ext rfl (snd_eq_one_of_mem_powTorsion hB x.2).symm)
  right_inv _ := rfl
  map_mul' _ _ := rfl

@[simp] theorem powTorsionProdEquiv_coe (m : ℕ) (hB : ∀ b : B, b ^ m = 1 → b = 1)
    (x : ↥(powTorsion (A × B) m)) :
    ((powTorsionProdEquiv (A := A) m hB x : ↥(powTorsion A m)) : A) = (x : A × B).1 := rfl

@[simp] theorem powTorsionProdEquiv_symm_coe (m : ℕ) (hB : ∀ b : B, b ^ m = 1 → b = 1)
    (a : ↥(powTorsion A m)) :
    (((powTorsionProdEquiv (A := A) (B := B) m hB).symm a : ↥(powTorsion (A × B) m)) : A × B)
      = ((a : A), 1) := rfl

/-- 捩れ部分群の元は `m` 乗すると `1`（部分型の中でも、もとの群の中でも）。 -/
theorem coe_pow_eq_one_of_mem_powTorsion {m : ℕ} (x : ↥(powTorsion A m)) : (x : A) ^ m = 1 := x.2

theorem pow_eq_one_of_mem_powTorsion {m : ℕ} (x : ↥(powTorsion A m)) : x ^ m = 1 :=
  Subtype.ext (by rw [SubmonoidClass.coe_pow, OneMemClass.coe_one]; exact x.2)

/-- ★**抽象核** —— 有限巡回群で全元が `m` 乗して `1` なら、位数は `m` を割る。

★これが「捩れの大きさは `p^n` を**割る**だけで、`p^n` とは限らない」の中身である。 -/
theorem natCard_dvd_of_isCyclic_of_forall_pow_eq_one {H : Type*} [Group H] [Finite H] [IsCyclic H]
    {m : ℕ} (h : ∀ x : H, x ^ m = 1) : Nat.card H ∣ m := by
  obtain ⟨g, hg⟩ := IsCyclic.exists_generator (α := H)
  have hcard : Nat.card H = orderOf g := (orderOf_eq_card_of_forall_mem_zpowers hg).symm
  rw [hcard]
  exact orderOf_dvd_of_pow_eq_one (h g)

end AbstractCore

/-! ## §2 `μ_m(K)` 側の一般論（体だけ。分岐は出てこない） -/

section RootsOfUnityCore

variable (F : Type*) [Field F]

/-- `μ_m(F)` の位数は `m` を割る。★`= m` ではない（`F` に依存する）。 -/
theorem natCard_rootsOfUnity_dvd (m : ℕ) [NeZero m] : Nat.card (rootsOfUnity m F) ∣ m := by
  refine natCard_dvd_of_isCyclic_of_forall_pow_eq_one (fun x => ?_)
  refine Subtype.ext (Units.ext ?_)
  have hx := (mem_rootsOfUnity m (x : Fˣ)).mp x.2
  push_cast
  simpa using congrArg (fun u : Fˣ => (u : F)) hx

end RootsOfUnityCore

/-! ## §3 具体層 1 —— `Gal(K^ab/K) ≅ 𝒪_K^× × Ẑ`

★Λ8 の `abelianGalEquivProd`（標的が `Gal(K^ur/K)`）を
`zhatMulEquivUnramGalArith` で `Ẑ` に取り替えるだけである。

★#59 回避: 定型 (b)。本ファイルには `restrictNormalHom` すら出てこない
（Λ8 が `Gal(K^ab/K)` の側で閉じた形にしてくれているので、
中間体の中の中間体を作る場面が 1 度も無い）。 -/

variable {p : ℕ} [Fact p.Prime]
variable (K : PAdicLocalField p)

section Decomposition

variable [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    {M : IntermediateField K.carrier K.closure}

/-- ★★**`Gal(M/K) ≅ 𝒪_K^× × Ẑ`**（`M = K_π · K^ur`）。

★Λ8 の `abelianGalEquivProd` と `zhatMulEquivUnramGalArith` を継ぐだけ（3 行）。 -/
noncomputable def galEquivUnitsZHat
    (hM : (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ⊔ unramifiedClosure K :
      IntermediateField K.carrier K.closure) = M) :
    (M ≃ₐ[K.carrier] M) ≃* ((𝒪[K.carrier])ˣ × ZHat) :=
  (abelianGalEquivProd K hq hπmax hπne0 f hf0 hf1 hf hM).trans
    (MulEquiv.prodCongr (MulEquiv.refl _) (zhatMulEquivUnramGalArith K).toMulEquiv.symm)

theorem galEquivUnitsZHat_apply
    (hM : (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ⊔ unramifiedClosure K :
      IntermediateField K.carrier K.closure) = M) (τ : M ≃ₐ[K.carrier] M) :
    galEquivUnitsZHat K hq hπmax hπne0 f hf0 hf1 hf hM τ
      = ((abelianGalEquivProd K hq hπmax hπne0 f hf0 hf1 hf hM τ).1,
          (zhatMulEquivUnramGalArith K).toMulEquiv.symm
            (abelianGalEquivProd K hq hπmax hπne0 f hf0 hf1 hf hM τ).2) := rfl

/-- ★★**位相群としての `Gal(M/K) ≃ₜ* 𝒪_K^× × Ẑ`**。

`continuous_abelianGalEquivProd` / `continuous_abelianGalEquivProd_symm`（Λ8）と
`ContinuousMulEquiv.prodCongr`（§2 の抽象核）を継ぐ。 -/
noncomputable def galContinuousEquivUnitsZHat [IsGalois K.carrier M]
    (hM : (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ⊔ unramifiedClosure K :
      IntermediateField K.carrier K.closure) = M) :
    (M ≃ₐ[K.carrier] M) ≃ₜ* ((𝒪[K.carrier])ˣ × ZHat) :=
  haveI := normal_lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf
  ContinuousMulEquiv.trans
    { toMulEquiv := abelianGalEquivProd K hq hπmax hπne0 f hf0 hf1 hf hM
      continuous_toFun := continuous_abelianGalEquivProd K hq hπmax hπne0 f hf0 hf1 hf hM
      continuous_invFun := continuous_abelianGalEquivProd_symm K hq hπmax hπne0 f hf0 hf1 hf hM }
    (ContinuousMulEquiv.prodCongr (ContinuousMulEquiv.refl _)
      (zhatMulEquivUnramGalArith K).symm)

theorem galContinuousEquivUnitsZHat_toMulEquiv [IsGalois K.carrier M]
    (hM : (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ⊔ unramifiedClosure K :
      IntermediateField K.carrier K.closure) = M) :
    (galContinuousEquivUnitsZHat K hq hπmax hπne0 f hf0 hf1 hf hM).toMulEquiv
      = galEquivUnitsZHat K hq hπmax hπne0 f hf0 hf1 hf hM := rfl

/-- `K^ab` への特殊化（★LKW を消費する段は Λ8 の
`abelianClosure_eq_lubinTateClosure_sup_unramifiedClosure` がやっている）。 -/
noncomputable def abelianGalEquivUnitsZHat :
    (abelianClosure K ≃ₐ[K.carrier] abelianClosure K) ≃* ((𝒪[K.carrier])ˣ × ZHat) :=
  galEquivUnitsZHat K hq hπmax hπne0 f hf0 hf1 hf
    (abelianClosure_eq_lubinTateClosure_sup_unramifiedClosure K hq hπmax hπne0 f hf0 hf1 hf).symm

/-- `K^ab` への特殊化（位相版）。 -/
noncomputable def abelianGalContinuousEquivUnitsZHat :
    (abelianClosure K ≃ₐ[K.carrier] abelianClosure K) ≃ₜ* ((𝒪[K.carrier])ˣ × ZHat) :=
  galContinuousEquivUnitsZHat K hq hπmax hπne0 f hf0 hf1 hf
    (abelianClosure_eq_lubinTateClosure_sup_unramifiedClosure K hq hπmax hπne0 f hf0 hf1 hf).symm

end Decomposition

/-- ★★**仮説なしの形** —— `Gal(K^ab/K) ≃ₜ* 𝒪_K^× × Ẑ`。

★同型は `π`（と Lubin-Tate 級数 `f`）の選択に依存するので `Nonempty` の内側に
閉じ込めてある。★**選択に依存しないとは主張していない。** -/
theorem nonempty_abelianGalContinuousEquivUnitsZHat (K : PAdicLocalField p) :
    Nonempty ((abelianClosure K ≃ₐ[K.carrier] abelianClosure K) ≃ₜ* ((𝒪[K.carrier])ˣ × ZHat)) := by
  haveI := isAdicComplete_valuationRing K
  haveI := valuationRing_isDVR K
  obtain ⟨ϖ, hϖirr⟩ := IsDiscreteValuationRing.exists_irreducible (𝒪[K.carrier])
  have hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {ϖ} :=
    (IsDiscreteValuationRing.irreducible_iff_uniformizer ϖ).mp hϖirr
  have hπne0 : ϖ ≠ 0 := hϖirr.ne_zero
  have hq : Fintype.card 𝓀[K.carrier] = p ^ (absoluteInertiaDegree K) := by
    rw [← Nat.card_eq_fintype_card]
    exact residueCard_eq_pow K
  obtain ⟨f, hf0, hf1, hf⟩ := exists_lubinTateSeries (A := 𝒪[K.carrier]) hq hπmax
  exact ⟨abelianGalContinuousEquivUnitsZHat K hq hπmax hπne0 f hf0 hf1 hf⟩

/-! ## §4 具体層 2 —— ★★到達点 `tors_m G ≅ μ_m(K)`

★`G` は「`𝒪_K^× × Ẑ` と同型な可換群」でありさえすればよい。
`G = Gal(K^ab/K)` はその 1 例である。 -/

section Torsion

variable {G : Type*} [CommGroup G] (e : G ≃* ((𝒪[K.carrier])ˣ × ZHat)) {m : ℕ} (hm : m ≠ 0)

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**(Λ9) `tors_m G ≅ μ_m(K)`**（`e : G ≅ 𝒪_K^× × Ẑ` を通して）。

在庫の `powTorsionCongr`（同型による捩れの移送）と
`torsionUnitsZHatEquiv`（`Ẑ` の捩れ自由性 + 「`m` 乗根は単数」）の合成。 -/
noncomputable def torsionEquivRootsOfUnity :
    ↥(powTorsion G m) ≃* ↥(rootsOfUnity m K.carrier) :=
  (powTorsionCongr e m).trans (torsionUnitsZHatEquiv K hm)

@[simp] theorem coe_torsionEquivRootsOfUnity (x : ↥(powTorsion G m)) :
    ((torsionEquivRootsOfUnity K e hm x : ↥(rootsOfUnity m K.carrier)) : (K.carrier)ˣ)
      = unitsToField K (e (x : G)).1 := rfl

/-- ★**`Ẑ` 因子は捩れに寄与しない** —— `tors_m G ≅ tors_m 𝒪_K^×`。

★§1 の抽象核 `powTorsionProdEquiv` に `Ẑ` の捩れ自由性
（`zhat_eq_one_of_pow_eq_one`）を代入しただけ。 -/
noncomputable def torsionEquivUnitsTorsion :
    ↥(powTorsion G m) ≃* ↥(powTorsion ((𝒪[K.carrier])ˣ) m) :=
  (powTorsionCongr e m).trans
    (powTorsionProdEquiv m (fun _ hb => zhat_eq_one_of_pow_eq_one hm hb))

include hm in
/-- 捩れ元の `Ẑ` 成分は自明。★#151 の実例——`hm` は結論に現れないので
`include` が要る（section の `variable` 仮説は theorem の本体だけからは拾われない）。
★`include … in` は docstring の**前**（#107/#154）。 -/
theorem snd_eq_one_of_mem_torsion (x : ↥(powTorsion G m)) : (e (x : G)).2 = 1 := by
  refine zhat_eq_one_of_pow_eq_one hm ?_
  have hx : (x : G) ^ m = 1 := x.2
  have h : e (x : G) ^ m = 1 := by rw [← map_pow, hx, map_one]
  exact congrArg Prod.snd h

end Torsion

/-! ## §5 作用 —— 大きさだけでは足りない

★`Gal(K^ab/K)` 自身の自己同型 `Ψ` が捩れの上で `u ↦ u^c` を誘導するなら、
`μ_m(K)` の側でも `ζ ↦ ζ^c` になる。★これが「作用が乗っている」の中身である。 -/

section Action

variable {G : Type*} [CommGroup G] (e : G ≃* ((𝒪[K.carrier])ˣ × ZHat)) {m : ℕ} (hm : m ≠ 0)

/-- ★**同型 `torsionEquivRootsOfUnity` は `c` 乗作用と可換**。

`Ψ` が `G` の（2 因子を混ぜてよい）自己同型で、捩れの上で第 1 成分が `c` 乗になるなら、
`μ_m(K)` の側でも `ζ ↦ ζ^c`。 -/
theorem torsionEquivRootsOfUnity_pow (Ψ : G ≃* G) (c : ℕ)
    (hΨ : ∀ x : G, x ^ m = 1 → (e (Ψ x)).1 = (e x).1 ^ c) (x : ↥(powTorsion G m)) :
    torsionEquivRootsOfUnity K e hm (powTorsionCongr Ψ m x)
      = (torsionEquivRootsOfUnity K e hm x) ^ c := by
  refine Subtype.ext ?_
  have hc : ((torsionEquivRootsOfUnity K e hm x ^ c : ↥(rootsOfUnity m K.carrier))
      : (K.carrier)ˣ) = (unitsToField K (e (x : G)).1) ^ c := by
    rw [SubmonoidClass.coe_pow, coe_torsionEquivRootsOfUnity]
  rw [coe_torsionEquivRootsOfUnity, hc, powTorsionCongr_coe, hΨ _ x.2, map_pow]

end Action

section ActionCyclotomic

variable {G : Type*} [CommGroup G]

/-- ★★**作用は円分指標倍**。

`K ⊆ F̄` を固定する。`g ∈ Γ_F` が `σ : K ≃ K` を誘導し、`G` の自己同型 `Ψ` が
捩れの上で `σ` を覆うなら、`tors_m G ≅ μ_m(K)` の下で `Ψ` は
`ζ ↦ ζ^{χ_{F,n}(g)}` として働く。指数の式は在庫の
`ProfiniteUnitsTorsion.lean::exists_torsionUnitsZHat_equiv_cyclotomic` と同じもの。 -/
theorem torsionEquivRootsOfUnity_cyclotomic (F : PAdicLocalField p) {n : ℕ}
    (ι : K.carrier →+* F.closure) (g : F.absGal) (σ : K.carrier ≃+* K.carrier)
    (hcompat : ∀ y : K.carrier, ι (σ y) = g (ι y))
    (e : G ≃* ((𝒪[K.carrier])ˣ × ZHat)) (Ψ : G ≃* G)
    (hΨ : ∀ x : G, x ^ (p ^ n) = 1 →
      ((unitsToField K (e (Ψ x)).1 : (K.carrier)ˣ) : K.carrier)
        = σ (((unitsToField K (e x).1 : (K.carrier)ˣ) : K.carrier)))
    (x : ↥(powTorsion G (p ^ n))) :
    torsionEquivRootsOfUnity K e (pn_ne_zero p n) (powTorsionCongr Ψ (p ^ n) x)
      = (torsionEquivRootsOfUnity K e (pn_ne_zero p n) x) ^
          ((PadicInt.toZModPow n
            ((cyclotomicCharacter F.closure p g.toRingEquiv : ℤ_[p]))).val) := by
  refine torsionEquivRootsOfUnity_pow K e (pn_ne_zero p n) Ψ _ (fun y hy => ?_) x
  have hu : (e y).1 ^ (p ^ n) = 1 := by
    have h : e y ^ (p ^ n) = 1 := by rw [← map_pow, hy, map_one]
    exact congrArg Prod.fst h
  refine units_eq_pow_of_coe K ?_
  rw [hΨ y hy, map_pow, Units.val_pow_eq_pow_val, unitsToField_coe,
    ringEquiv_pow_eq_cyclotomicCharacter F K ι g σ hcompat
      (coe_pow_eq_one_of_units_pow_eq_one K hu)]

end ActionCyclotomic

/-! ## §6 ★★到達点 —— `tors_{p^n}(Gal(K^ab/K)) ≅ μ_{p^n}(K)` -/

section AbelianGal

variable [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))

/-- **`tors_{p^n}(Gal(K^ab/K)) ≅ μ_{p^n}(K)`**（素元・級数を明示した形）。 -/
noncomputable def abelianGalTorsionEquivRootsOfUnity (n : ℕ) :
    ↥(powTorsion (abelianClosure K ≃ₐ[K.carrier] abelianClosure K) (p ^ n))
      ≃* ↥(rootsOfUnity (p ^ n) K.carrier) :=
  torsionEquivRootsOfUnity K (abelianGalEquivUnitsZHat K hq hπmax hπne0 f hf0 hf1 hf)
    (pn_ne_zero p n)

end AbelianGal

/-- 捩れ部分群そのものは Lubin-Tate の選択に依存しない
（`{x | x^{p^n} = 1}` そのものである）。★依存するのは同型の方だけ。 -/
theorem mem_powTorsion_abelianGal (K : PAdicLocalField p) {n : ℕ}
    {x : abelianClosure K ≃ₐ[K.carrier] abelianClosure K} :
    x ∈ powTorsion (abelianClosure K ≃ₐ[K.carrier] abelianClosure K) (p ^ n)
      ↔ x ^ (p ^ n) = 1 := Iff.rfl

def exists_abelianGalTorsion_equiv_rootsOfUnity.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**(Λ9、仮説ゼロ) `tors_{p^n}(Gal(K^ab/K)) ≅ μ_{p^n}(K)`。**

原文 (pGC p.3):
> Now recall from local class field theory (see, e.g., [3]) that we have a natural isomorphism ΓKab ≅ (K×)∧ where the superscripted "ab" denotes the abelianization, and "(K×)∧" denotes the profinite completion of K× = K − {0}. Let k be the residue field of OK (the ring of integers of K). Thus, k is the field of q = pf elements. Now it is well-known that (K×)∧ fits into an exact sequence of topological groups: 0 → UK → (K×)∧ → Z → 0 where UK d=ef OK×.

★★これが pGC Proposition 1.1 が名指しで推奨する道の最後の環である。

★同型 `e` は `∃`（`Nonempty`）の内側にあり、結論の型は `K` と `n` にしか依存しない。
★★**右辺は `K` に依存する群**であって `ZMod (p^n)` ではない
（`K = ℚ_p`、`p` 奇 なら `μ_p(ℚ_p) = 1` で自明群になる）。 -/
theorem exists_abelianGalTorsion_equiv_rootsOfUnity (K : PAdicLocalField p) (n : ℕ) :
    Nonempty (↥(powTorsion (abelianClosure K ≃ₐ[K.carrier] abelianClosure K) (p ^ n))
      ≃* ↥(rootsOfUnity (p ^ n) K.carrier)) := by
  obtain ⟨Φ⟩ := nonempty_abelianGalContinuousEquivUnitsZHat K
  exact ⟨torsionEquivRootsOfUnity K Φ.toMulEquiv (pn_ne_zero p n)⟩

/-- ★**`Ẑ` 因子は捩れに寄与しない**（`K^ab` 版）——
`tors_{p^n}(Gal(K^ab/K)) ≅ tors_{p^n}(𝒪_K^×)`。 -/
theorem exists_abelianGalTorsion_equiv_unitsTorsion (K : PAdicLocalField p) (n : ℕ) :
    Nonempty (↥(powTorsion (abelianClosure K ≃ₐ[K.carrier] abelianClosure K) (p ^ n))
      ≃* ↥(powTorsion ((𝒪[K.carrier])ˣ) (p ^ n))) := by
  obtain ⟨Φ⟩ := nonempty_abelianGalContinuousEquivUnitsZHat K
  exact ⟨torsionEquivUnitsTorsion K Φ.toMulEquiv (pn_ne_zero p n)⟩

def exists_abelianGalTorsion_equiv_cyclotomic.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**(Λ9、作用つき) `tors_{p^n}(Gal(K^ab/K)) ≅ μ_{p^n}(K)` は円分指標倍を運ぶ。**

`K ⊆ F̄` の中で `g ∈ Γ_F` が `σ : K ≃ K` を誘導し、`Gal(K^ab/K)` の自己同型 `Ψ` が
捩れの上で `σ` を覆うなら、`e` の下で `Ψ` は `ζ ↦ ζ^{χ_{F,n}(g)}` として働く。

★大きさだけでは `TorsionCyclotomeIsCyclotomic` に接続できない——
接続に要るのはこの**指数の式**である。 -/
theorem exists_abelianGalTorsion_equiv_cyclotomic (F K : PAdicLocalField p) (n : ℕ)
    (ι : K.carrier →+* F.closure) :
    ∃ (E : (abelianClosure K ≃ₐ[K.carrier] abelianClosure K) ≃* ((𝒪[K.carrier])ˣ × ZHat))
      (e : ↥(powTorsion (abelianClosure K ≃ₐ[K.carrier] abelianClosure K) (p ^ n))
        ≃* ↥(rootsOfUnity (p ^ n) K.carrier)),
      (∀ x, ((e x : ↥(rootsOfUnity (p ^ n) K.carrier)) : (K.carrier)ˣ)
          = unitsToField K (E (x : abelianClosure K ≃ₐ[K.carrier] abelianClosure K)).1)
      ∧ ∀ (g : F.absGal) (σ : K.carrier ≃+* K.carrier)
          (Ψ : (abelianClosure K ≃ₐ[K.carrier] abelianClosure K)
            ≃* (abelianClosure K ≃ₐ[K.carrier] abelianClosure K)),
          (∀ y : K.carrier, ι (σ y) = g (ι y)) →
          (∀ x : (abelianClosure K ≃ₐ[K.carrier] abelianClosure K), x ^ (p ^ n) = 1 →
            ((unitsToField K (E (Ψ x)).1 : (K.carrier)ˣ) : K.carrier)
              = σ (((unitsToField K (E x).1 : (K.carrier)ˣ) : K.carrier))) →
          ∀ x, e (powTorsionCongr Ψ (p ^ n) x)
            = (e x) ^ ((PadicInt.toZModPow n
                ((cyclotomicCharacter F.closure p g.toRingEquiv : ℤ_[p]))).val) := by
  obtain ⟨Φ⟩ := nonempty_abelianGalContinuousEquivUnitsZHat K
  exact ⟨Φ.toMulEquiv, torsionEquivRootsOfUnity K Φ.toMulEquiv (pn_ne_zero p n),
    fun _ => rfl,
    fun g σ Ψ hcompat hΨ x =>
      torsionEquivRootsOfUnity_cyclotomic K F ι g σ hcompat Φ.toMulEquiv Ψ hΨ x⟩

/-! ## §7 退化の自己検査 -/

section Degeneracy

/-- `tors_{p^n}(Gal(K^ab/K))` は**有限**。 -/
theorem finite_powTorsion_abelianGal (K : PAdicLocalField p) (n : ℕ) :
    Finite ↥(powTorsion (abelianClosure K ≃ₐ[K.carrier] abelianClosure K) (p ^ n)) := by
  haveI : NeZero (p ^ n) := ⟨pn_ne_zero p n⟩
  obtain ⟨e⟩ := exists_abelianGalTorsion_equiv_rootsOfUnity K n
  exact Finite.of_equiv _ e.toEquiv.symm

/-- `tors_{p^n}(Gal(K^ab/K))` は**巡回**。★`μ_{p^n}(K)` が体の中の 1 の冪根の群だから。 -/
theorem isCyclic_powTorsion_abelianGal (K : PAdicLocalField p) (n : ℕ) :
    IsCyclic ↥(powTorsion (abelianClosure K ≃ₐ[K.carrier] abelianClosure K) (p ^ n)) := by
  haveI : NeZero (p ^ n) := ⟨pn_ne_zero p n⟩
  obtain ⟨e⟩ := exists_abelianGalTorsion_equiv_rootsOfUnity K n
  exact (MulEquiv.isCyclic e).mpr (rootsOfUnity.isCyclic K.carrier (p ^ n))

/-- ★★**位数は `p^n` を割る。★`= p^n` とは限らない。**

`K = ℚ_p`（`p` 奇）なら `μ_p(ℚ_p) = 1` なので `n = 1` で位数 `1` になる。
★「常に位数 `p^n`」と書くのは誤りである。 -/
theorem natCard_powTorsion_abelianGal_dvd (K : PAdicLocalField p) (n : ℕ) :
    Nat.card ↥(powTorsion (abelianClosure K ≃ₐ[K.carrier] abelianClosure K) (p ^ n)) ∣ p ^ n := by
  haveI : NeZero (p ^ n) := ⟨pn_ne_zero p n⟩
  haveI := finite_powTorsion_abelianGal K n
  haveI := isCyclic_powTorsion_abelianGal K n
  refine natCard_dvd_of_isCyclic_of_forall_pow_eq_one (fun x => ?_)
  exact pow_eq_one_of_mem_powTorsion x

/-- ★**`μ_{p^n}(K)` の位数も `p^n` を割るだけ**（同じ理由）。 -/
theorem natCard_rootsOfUnity_pow_dvd (K : PAdicLocalField p) (n : ℕ) :
    Nat.card (rootsOfUnity (p ^ n) K.carrier) ∣ p ^ n := by
  haveI : NeZero (p ^ n) := ⟨pn_ne_zero p n⟩
  exact natCard_rootsOfUnity_dvd K.carrier (p ^ n)

end Degeneracy

end ABC3.Found.PGC
