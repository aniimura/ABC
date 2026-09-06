import ABC3.Found.PGC.AbelianSplitUnramified

/-!
# Milne が黙っている段——Frobenius の持ち上げを `Gal(K̄/F)` の中から選ぶ

Milne, *Class Field Theory*, p.57 の THEOREM 4.8 の証明は、Lemma 4.11 が作る
完全分岐部分体 `L_t` について「(4.9) implies that `L_t ⊆ K_π`」と書く。
しかし Lemma 4.9 の証明は

> PROOF. Let `G = Gal(L/K)` and `H = Gal(L/K_π)`, so that `G/H = Gal(K_π/K)`.

と始まる。`H = Gal(L/K_π)` が `G` の部分群であり `G/H = Gal(K_π/K)` であるためには
**`K_π ⊆ L`** でなければならない。すなわち Lemma 4.9 の仮定は「`L ⊇ K_π`」であり、
4.8 の証明でこれを `L_t` に当てるには **`K_π ⊆ L_t`** が要る
(原文の "If `L ⊃ K_π`, then `L = K_π`" の `⊃` は pdftotext で消える文字だが、
証明の 1 行目が向きを確定させる)。

★**この包含は原文からは出てこない。** Lemma 4.11 の証明は
「`σ ∈ Gal(LK_m/K)` を `σ|_{K_m}` が Frobenius になるように取る」としか言わず、
`σ` が `K_π` を固定するとは指定していない。`L_t = (LK_m)^{⟨σ⟩}` が `K_π` を含むのは
**`σ` を `Gal(K̄/K_π)` の中から選んだとき**に限る。選べる理由は
`K_π ⊓ K_m = ⊥`(完全分岐 ⊓ 不分岐)——このとき制限写像
`Gal(K̄/K_π) → Gal(K_m/K)` が全射になる。

本ファイルはその選択を明示的に行い、`Found/PGC/AbelianFrobeniusSplit.lean` と
`Found/PGC/AbelianSplitUnramified.lean` の主定理の**強化版**を与える。
★既存ファイルは一切書き換えていない(import のみ)。

## 主定理

* `exists_unramified_frobenius_lift_fixing`——`F ⊓ K^ur = ⊥` なる中間体 `F` に対し、
  `exists_unramified_frobenius_lift` の `σ` を **`σ ∈ F.fixingSubgroup`**
  (= `σ` は `F` を各点固定する)に取れる。
* `exists_totallyRamified_abelian_split_ge`——Milne Lemma 4.11 の強化版。
  `F ≤ L` かつ `F ⊓ K^ur = ⊥` なら、`L_t = K(α)` を **`F ≤ K(α)`** に取れる。
* `exists_totallyRamified_abelian_split_ge_unramLevel`——同じものを
  `unramLevel K m` で名指しした形。

## 全射性に使った道具(在庫にあった)

`AbelianDecomposition.lean::sup_fixingSubgroup_eq_top`
(`A ⊓ B = ⊥`・`B/k` 正規 ⟹ `Gal(E/A) ⊔ Gal(E/B) = ⊤`)が**そのまま**使えた。
これと `Subgroup.mul_normal`(`↑(A ⊔ N) = ↑A * ↑N`)を合わせると
`Γ_K = Gal(K̄/F) · Gal(K̄/K_m)` となり、任意の持ち上げ `σ₀ = σ · n`
(`σ ∈ Gal(K̄/F)`,`n ∈ Gal(K̄/K_m)`)から `F` を固定する持ち上げ `σ` が取れる。
`σ` と `σ₀` は `Gal(K̄/K_m)` を法として等しいので、`AbelianFrobeniusSplit.lean` の
`σ` に課された 4 条件はすべて `σ` へ移る(`zpow_mem_iff_of_inv_mul_mem` /
`zpowers_sup_eq_top_of_inv_mul_mem`——どちらも純群論、商群 `Γ/N` を 1 回だけ使う)。

## ★退化の自己検査

* **`F ⊓ K^ur = ⊥` を落とすと偽。** 具体例:`K = ℚ_p`、`m = 2`、
  `F = K_2`(次数 2 の不分岐拡大、これは `K_m` そのもの)。
  `σ ∈ F.fixingSubgroup` は `K_2` を各点固定するので `Gal(K_2/K) ≅ ℤ/2` への像は
  単位元であり、生成元(Frobenius)には決してならない。すなわち
  「`F` を固定する Frobenius 持ち上げ」は存在しない。
  一般に `F ⊓ K_m ≠ ⊥` なら `Gal(K̄/F) → Gal(K_m/K)` の像は
  `Gal(K_m/(F ⊓ K_m))` に含まれ、真部分群になる。
* **`exists_totallyRamified_abelian_split_ge` でも同じ仮定は落とせない。**
  `F ≤ K(α)` と `K(α)/K` 完全分岐から `F ⊓ K^ur ≤ K(α) ⊓ K^ur = ⊥` が
  従うので、結論が仮定を含意している(仮定は必要十分の側)。
* **`F ≤ L` も落とせない。** `F ⊄ L ⊔ K_m` なら `L_t ⊔ K_m = L ⊔ K_m` と
  `F ≤ L_t` は両立しない。
* ★結論に自由なパラメータを出していない(`w`・`σ`・`α`・`β` はすべて `∃` の内側)。

## 逸脱(記録)

1. **仮定を `F ⊓ K_m = ⊥` でなく `F ⊓ K^ur = ⊥` にした**(強めた)。
   `K_m` は結論の `∃` の内側で作られるので仮定に書けない。
   `K_m ≤ K^ur` なので `F ⊓ K^ur = ⊥ ⟹ F ⊓ K_m = ⊥`。
   `F` が `K` 上完全分岐なら `totallyRamifiedAdjoin_inf_unramifiedClosure` が
   この仮定をそのまま供給するので、下流(`F = K_{π,n}`)では負担にならない。
2. **`F` は有限次に限られる**(`F ≤ L` かつ `L` は有限次)。原文が
   Lemma 4.11 を当てている `L · K_π` は無限次で、Milne の `L_t` も
   無限次(`L_t = K_π`)になる。本ファイルは**有限段** `F = K_{π,n} ≤ L · K_{π,n}`
   だけを扱う。無限次への極限(`K_π = ⋃ K_{π,n}` と `L_t` の増大列)は
   本ファイルの担当ではない。★したがって「4.8 の証明の穴を塞ぐ」うちの
   **持ち上げの選択**は塞いだが、**無限段への移行**は残っている。
3. 原文の `σ|_{K_m} = Frobenius` は、`AbelianFrobeniusSplit.lean` と同じく
   「`Gal(K_m/K)` の生成元の持ち上げ」として実現している(剰余体の `x ↦ x^q`
   との一致は主張しない)。
4. `exists_totallyRamified_abelian_split_ge` の証明は
   `AbelianSplitUnramified.lean::exists_totallyRamified_abelian_split` の証明を
   なぞっている(既存ファイルを書き換えない、という制約のため複製した)。
   差分は 2 か所だけ:`σ` を `exists_unramified_frobenius_lift_fixing` から取ることと、
   最後に `S ≤ Gal(K̄/F)` から `F ≤ L_t` を出すこと。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued
open Pointwise

/-! ## 0. 「正規部分群を法として等しい 2 元」の乗り換え(純群論) -/

/-- `A ⊔ N = ⊤`(`N ◁ Γ`)なら、任意の `x` は `a * n`(`a ∈ A`,`n ∈ N`)と書ける。

`Subgroup.mul_normal`(`↑(A ⊔ N) = ↑A * ↑N`)を `⊤` に当てるだけ。 -/
theorem exists_mem_mul_mem_of_sup_eq_top {Γ : Type*} [Group Γ] {A N : Subgroup Γ} [N.Normal]
    (h : A ⊔ N = ⊤) (x : Γ) : ∃ a ∈ A, ∃ n ∈ N, x = a * n := by
  have hx : x ∈ ((A : Set Γ) * (N : Set Γ)) := by
    rw [← Subgroup.mul_normal, h]; trivial
  obtain ⟨a, ha, n, hn, hx⟩ := hx
  exact ⟨a, ha, n, hn, hx.symm⟩

/-- `σ` と `h` が正規部分群 `N` を法として等しければ、**冪の `N` 所属は一致する**。

商群 `Γ ⧸ N` を 1 回だけ作って `map_zpow` で移す。 -/
theorem zpow_mem_iff_of_inv_mul_mem {Γ : Type*} [Group Γ] {N : Subgroup Γ} [N.Normal] {σ h : Γ}
    (hmem : σ⁻¹ * h ∈ N) (k : ℤ) : h ^ k ∈ N ↔ σ ^ k ∈ N := by
  have hq : (QuotientGroup.mk' N) σ = (QuotientGroup.mk' N) h := by
    simp only [QuotientGroup.mk'_apply]
    exact QuotientGroup.eq.mpr hmem
  constructor
  · intro hk
    have hh : (QuotientGroup.mk' N) (h ^ k) = 1 := (QuotientGroup.eq_one_iff _).mpr hk
    rw [map_zpow, ← hq, ← map_zpow] at hh
    exact (QuotientGroup.eq_one_iff _).mp hh
  · intro hk
    have hh : (QuotientGroup.mk' N) (σ ^ k) = 1 := (QuotientGroup.eq_one_iff _).mpr hk
    rw [map_zpow, hq, ← map_zpow] at hh
    exact (QuotientGroup.eq_one_iff _).mp hh

/-- `N` を法として `σ` と等しい `h` について `⟨σ⟩ ⊔ N = ⊤ ⟹ ⟨h⟩ ⊔ N = ⊤`。

`σ = h * (σ⁻¹h)⁻¹ ∈ ⟨h⟩ ⊔ N` なので `⟨σ⟩ ⊔ N ≤ ⟨h⟩ ⊔ N`。
★ここは `N` の正規性を使わない。 -/
theorem zpowers_sup_eq_top_of_inv_mul_mem {Γ : Type*} [Group Γ] {N : Subgroup Γ} {σ h : Γ}
    (hmem : σ⁻¹ * h ∈ N) (htop : Subgroup.zpowers σ ⊔ N = ⊤) :
    Subgroup.zpowers h ⊔ N = ⊤ := by
  refine top_le_iff.mp ?_
  rw [← htop]
  refine sup_le ?_ le_sup_right
  rw [Subgroup.zpowers_le]
  have hσ : σ = h * (σ⁻¹ * h)⁻¹ := by group
  rw [hσ]
  exact Subgroup.mul_mem _
    (le_sup_left (a := Subgroup.zpowers h) (b := N) (Subgroup.mem_zpowers h))
    (le_sup_right (a := Subgroup.zpowers h) (b := N) (Subgroup.inv_mem N hmem))

variable {p : ℕ} [Fact p.Prime]

/-! ## 1. `F` を固定する Frobenius 持ち上げ -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**Milne が指定していない選択**——`F ⊓ K^ur = ⊥` なる中間体 `F` に対し、
`exists_unramified_frobenius_lift` の Frobenius 持ち上げ `σ` を
**`F` を各点固定するもの**に取れる。

すなわち次数 `m` の不分岐拡大 `K_m = K(w)` と `σ ∈ Γ_K` が在って

* `σ ∈ Gal(K̄/F)`(★これが原文に無い保証)
* `σ^k ∈ Gal(K̄/K_m) ⟺ m ∣ k`
* `⟨σ⟩ ⊔ Gal(K̄/K_m) = ⊤`
* `∀ g, g^m ∈ Gal(K̄/K_m)`、`∀ a b, [a,b] ∈ Gal(K̄/K_m)`

を満たす。

## 証明

`F ⊓ K_m = ⊥`(`K_m ≤ K^ur` と仮定から)なので
`sup_fixingSubgroup_eq_top` により `Gal(K̄/F) ⊔ Gal(K̄/K_m) = ⊤`。
`Gal(K̄/K_m)` は正規(`normal_fixingSubgroup`)だから
`Γ_K = Gal(K̄/F) · Gal(K̄/K_m)`。任意の持ち上げ `σ₀` を `σ · n` と分解すれば
`σ ∈ Gal(K̄/F)` で、`σ₀⁻¹σ = n⁻¹ ∈ Gal(K̄/K_m)` ゆえ残りの条件は
`zpow_mem_iff_of_inv_mul_mem` / `zpowers_sup_eq_top_of_inv_mul_mem` で移る。

退化の自己検査:`F ⊓ K^ur = ⊥` を落とすと**偽**。`K = ℚ_p`、`m = 2`、
`F = K_2`(次数 2 の不分岐拡大)とすると、`F` を固定する `σ` の
`Gal(K_2/K) ≅ ℤ/2` への像は単位元で、生成元にならない
(`σ^1 ∈ Gal(K̄/K_2)` かつ `¬ (2 ∣ 1)` で第 2 条件が破れる)。 -/
theorem exists_unramified_frobenius_lift_fixing (K : PAdicLocalField p) {m : ℕ} (hm : m ≠ 0)
    (F : IntermediateField K.carrier K.closure)
    (hF : F ⊓ unramifiedClosure K = ⊥) :
    ∃ (w : K.closure) (σ : K.absGal),
      IsUnramifiedAdjoin K w ∧
      Module.finrank K.carrier
        (IntermediateField.adjoin K.carrier ({w} : Set K.closure)) = m ∧
      σ ∈ F.fixingSubgroup ∧
      (∀ k : ℤ, σ ^ k ∈ (IntermediateField.adjoin K.carrier
          ({w} : Set K.closure)).fixingSubgroup ↔ (m : ℤ) ∣ k) ∧
      Subgroup.zpowers σ ⊔ (IntermediateField.adjoin K.carrier
        ({w} : Set K.closure)).fixingSubgroup = ⊤ ∧
      (∀ g : K.absGal, g ^ m ∈ (IntermediateField.adjoin K.carrier
        ({w} : Set K.closure)).fixingSubgroup) ∧
      (∀ a b : K.absGal, a * b * a⁻¹ * b⁻¹ ∈ (IntermediateField.adjoin K.carrier
        ({w} : Set K.closure)).fixingSubgroup) := by
  haveI := isGalois_closure K
  obtain ⟨w, σ0, huw, hrankw, hkm, htopσ, hpowm, hcommm⟩ := exists_unramified_frobenius_lift K hm
  haveI := normal_of_isUnramifiedAdjoin K w huw
  set Km := IntermediateField.adjoin K.carrier ({w} : Set K.closure) with hKmdef
  haveI : (Km.fixingSubgroup : Subgroup K.absGal).Normal := normal_fixingSubgroup Km
  have hFKm : F ⊓ Km = ⊥ :=
    le_antisymm (le_trans (inf_le_inf_left F (adjoin_le_unramifiedClosure K huw)) hF.le) bot_le
  have htopF : F.fixingSubgroup ⊔ Km.fixingSubgroup = ⊤ :=
    sup_fixingSubgroup_eq_top F Km hFKm
  obtain ⟨σ, hσF, n, hn, hσ0⟩ := exists_mem_mul_mem_of_sup_eq_top htopF σ0
  have hmem : σ0⁻¹ * σ ∈ (Km.fixingSubgroup : Subgroup K.absGal) := by
    have hrw : σ0⁻¹ * σ = n⁻¹ := by rw [hσ0]; group
    rw [hrw]
    exact Subgroup.inv_mem _ hn
  refine ⟨w, σ, huw, hrankw, hσF, ?_, ?_, hpowm, hcommm⟩
  · intro k
    rw [zpow_mem_iff_of_inv_mul_mem hmem k]
    exact hkm k
  · exact zpowers_sup_eq_top_of_inv_mul_mem hmem htopσ

/-! ## 2. Milne Lemma 4.11 の強化版(`F ≤ L_t`) -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**Milne Lemma 4.11(`L_t` が `F` を含む版)**——
`L/K` が有限アーベルで `Monoid.exponent Gal(L/K) ∣ m`(`m ≠ 0`)、
さらに `F ≤ L` かつ `F ⊓ K^ur = ⊥` のとき、`K.closure` の中に

* `α`:`K(α)/K` は**完全分岐**、`K` 上正規、`Gal(K(α)/K)` は可換、
  ★しかも **`F ≤ K(α)`**
* `β`:`K(β)/K` は**不分岐で次数 `m`**(すなわち `K_m`)

が在って `K(α) ⊔ K(β) = L ⊔ K(β)`。

`F ≤ K(α)` が Milne の THEOREM 4.8 の証明で Lemma 4.9 を当てるのに要る包含
(`K_π ⊆ L_t`)の有限段版である。

証明は `exists_totallyRamified_abelian_split` と同じ 4 段で、`σ` を
`exists_unramified_frobenius_lift_fixing` から取る点だけが違う。すると
`S = ⟨σ⟩ ⊔ Gal(K̄/(L·K_m)) ≤ Gal(K̄/F)`(第 1 項は `σ ∈ Gal(K̄/F)`、
第 2 項は `F ≤ L ≤ L·K_m`)なので `F = fixedField Gal(K̄/F) ≤ fixedField S = L_t`。

★退化:`F ⊓ K^ur = ⊥` は落とせない(`F ≤ K(α)` と `K(α)/K` 完全分岐から
逆に従うので、仮定と結論が同値な部分)。`F ≤ L` も落とせない。 -/
theorem exists_totallyRamified_abelian_split_ge (K : PAdicLocalField p)
    (L : IntermediateField K.carrier K.closure) [FiniteDimensional K.carrier L]
    [Normal K.carrier L] (hab : ∀ a b : (L ≃ₐ[K.carrier] L), a * b = b * a)
    {m : ℕ} (hm : m ≠ 0)
    (hexp : Monoid.exponent (L ≃ₐ[K.carrier] L) ∣ m)
    (F : IntermediateField K.carrier K.closure) (hFL : F ≤ L)
    (hFur : F ⊓ unramifiedClosure K = ⊥) :
    ∃ α β : K.closure,
      IsTotallyRamifiedAdjoin K α ∧
      Normal K.carrier (IntermediateField.adjoin K.carrier ({α} : Set K.closure)) ∧
      (∀ a b : ((IntermediateField.adjoin K.carrier ({α} : Set K.closure))
        ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({α} : Set K.closure))),
          a * b = b * a) ∧
      F ≤ IntermediateField.adjoin K.carrier ({α} : Set K.closure) ∧
      IsUnramifiedAdjoin K β ∧
      Module.finrank K.carrier
        (IntermediateField.adjoin K.carrier ({β} : Set K.closure)) = m ∧
      IntermediateField.adjoin K.carrier ({α} : Set K.closure)
          ⊔ IntermediateField.adjoin K.carrier ({β} : Set K.closure)
        = L ⊔ IntermediateField.adjoin K.carrier ({β} : Set K.closure) := by
  haveI := isGalois_closure K
  haveI := IsAlgClosure.normal K.carrier K.closure
  obtain ⟨w, σ, huw, hrankw, hσF, hkm, htopσ, hpowm, hcommm⟩ :=
    exists_unramified_frobenius_lift_fixing K hm F hFur
  haveI := normal_of_isUnramifiedAdjoin K w huw
  set Km := IntermediateField.adjoin K.carrier ({w} : Set K.closure) with hKmdef
  haveI : (L.fixingSubgroup : Subgroup K.absGal).Normal := normal_fixingSubgroup L
  haveI : (Km.fixingSubgroup : Subgroup K.absGal).Normal := normal_fixingSubgroup Km
  have hMfix : (L ⊔ Km).fixingSubgroup = L.fixingSubgroup ⊓ Km.fixingSubgroup :=
    fixingSubgroup_sup L Km
  -- 第 2 段:`Γ_K/H_M` の指数は `m` を割る
  have hexpM : ∀ g : K.absGal, g ^ m ∈ L.fixingSubgroup ⊓ Km.fixingSubgroup := by
    intro g
    have h1 : g ^ m ∈ (L.fixingSubgroup : Subgroup K.absGal) := by
      rw [← ker_restrictNormalHom_eq_fixingSubgroup L, MonoidHom.mem_ker, map_pow]
      obtain ⟨c, hc⟩ := hexp
      rw [hc, pow_mul, Monoid.pow_exponent_eq_one, one_pow]
    exact Subgroup.mem_inf.mpr ⟨h1, hpowm g⟩
  have hcommM : ∀ a b : K.absGal,
      a * b * a⁻¹ * b⁻¹ ∈ L.fixingSubgroup ⊓ Km.fixingSubgroup := by
    intro a b
    refine Subgroup.mem_inf.mpr ⟨?_, hcommm a b⟩
    have h1 : (AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure) (L : Type _))
        (a * b * a⁻¹ * b⁻¹) = 1 := by
      simp only [map_mul, map_inv]
      exact commutator_eq_one_of_mul_comm (hab _ _)
    rw [← ker_restrictNormalHom_eq_fixingSubgroup L, MonoidHom.mem_ker]
    exact h1
  have hHMm : (L.fixingSubgroup ⊓ Km.fixingSubgroup : Subgroup K.absGal) ≤ Km.fixingSubgroup :=
    inf_le_right
  have hHMopen : IsOpen ((L.fixingSubgroup ⊓ Km.fixingSubgroup : Subgroup K.absGal)
      : Set K.absGal) := by
    rw [Subgroup.coe_inf]
    exact (IntermediateField.fixingSubgroup_isOpen L).inter
      (IntermediateField.fixingSubgroup_isOpen Km)
  -- 第 3 段:`S := ⟨σ⟩ ⊔ H_M` と、その固定体 `L_t`
  set S := Subgroup.zpowers σ ⊔ (L.fixingSubgroup ⊓ Km.fixingSubgroup) with hSdef
  have hSclosed : IsClosed (S : Set K.absGal) :=
    Subgroup.isClosed_of_isOpen _ (Subgroup.isOpen_mono le_sup_right hHMopen)
  set Lt := IntermediateField.fixedField S with hLtdef
  have hfixLt : Lt.fixingSubgroup = S :=
    InfiniteGalois.fixingSubgroup_fixedField (⟨S, hSclosed⟩ : ClosedSubgroup K.absGal)
  have hσk : ∀ k : ℤ, σ ^ k ∈ (Km.fixingSubgroup : Subgroup K.absGal) →
      σ ^ k ∈ (L.fixingSubgroup ⊓ Km.fixingSubgroup : Subgroup K.absGal) := by
    intro k hk
    obtain ⟨c, hc⟩ := (hkm k).mp hk
    rw [hc, zpow_mul, zpow_natCast]
    exact Subgroup.zpow_mem _ (hexpM σ) c
  have hSHm : S ⊓ Km.fixingSubgroup = L.fixingSubgroup ⊓ Km.fixingSubgroup :=
    zpowers_sup_inf_eq hHMm hσk
  have hsup : Lt ⊔ Km = L ⊔ Km := by
    have h1 : (Lt ⊔ Km).fixingSubgroup = (L ⊔ Km).fixingSubgroup := by
      rw [fixingSubgroup_sup, hfixLt, hSHm, hMfix]
    have h2 := congrArg IntermediateField.fixedField h1
    rwa [InfiniteGalois.fixedField_fixingSubgroup,
      InfiniteGalois.fixedField_fixingSubgroup] at h2
  have hfixedTop : IntermediateField.fixedField (⊤ : Subgroup K.absGal)
      = (⊥ : IntermediateField K.carrier K.closure) := by
    rw [← IntermediateField.fixingSubgroup_bot (F := K.carrier) (E := K.closure)]
    exact InfiniteGalois.fixedField_fixingSubgroup ⊥
  have hStop : S ⊔ Km.fixingSubgroup = ⊤ := by
    rw [hSdef, sup_assoc, sup_eq_right.mpr hHMm]
    exact htopσ
  have hinf : Lt ⊓ Km = ⊥ := by
    have h := fixedField_sup S (Km.fixingSubgroup)
    rw [hStop, hfixedTop, InfiniteGalois.fixedField_fixingSubgroup Km] at h
    exact h.symm
  -- 第 4 段:`L ⊔ K_m` の不分岐部分は `K_m` に入る
  have hMur : ((L ⊔ Km) ⊓ unramifiedClosure K : IntermediateField K.carrier K.closure) ≤ Km := by
    intro z hz
    obtain ⟨hzM, hzur⟩ := hz
    have huz : IsUnramifiedAdjoin K z := (mem_unramifiedClosure_iff_isUnramified K z).mp hzur
    have hzle : IntermediateField.adjoin K.carrier ({z} : Set K.closure) ≤ L ⊔ Km := by
      rw [IntermediateField.adjoin_simple_le_iff]; exact hzM
    have hPle : (L.fixingSubgroup ⊓ Km.fixingSubgroup : Subgroup K.absGal)
        ≤ (IntermediateField.adjoin K.carrier ({z} : Set K.closure)).fixingSubgroup := by
      rw [← hMfix]; exact IntermediateField.fixingSubgroup_le hzle
    have hdvd := finrank_dvd_of_isUnramified_of_pow_mem K huz (fun g => hPle (hexpM g))
    exact adjoin_le_of_dvd K z w huz huw (by rw [hrankw]; exact hdvd)
      (IntermediateField.mem_adjoin_simple_self _ _)
  have hLtle : Lt ≤ L ⊔ Km := by rw [← hsup]; exact le_sup_left
  have hLtur : Lt ⊓ unramifiedClosure K = ⊥ := by
    refine le_antisymm ?_ bot_le
    refine le_trans (le_inf inf_le_left ?_) hinf.le
    exact le_trans (inf_le_inf hLtle le_rfl) hMur
  -- `L_t/K` のアーベル性(ここだけが `hab` を使う)
  have hSnormal : S.Normal := by
    constructor
    intro s hs g
    have h1 : g * s * g⁻¹ * s⁻¹ ∈ (L.fixingSubgroup ⊓ Km.fixingSubgroup : Subgroup K.absGal) :=
      hcommM g s
    have h2 : g * s * g⁻¹ = (g * s * g⁻¹ * s⁻¹) * s := by group
    rw [h2]
    exact Subgroup.mul_mem _ (le_sup_right (a := Subgroup.zpowers σ) h1) hs
  haveI : (Lt.fixingSubgroup : Subgroup K.absGal).Normal := by rw [hfixLt]; exact hSnormal
  haveI : IsGalois K.carrier Lt := (InfiniteGalois.normal_iff_isGalois Lt).mp inferInstance
  have hLtcomm : ∀ a b : (Lt ≃ₐ[K.carrier] Lt), a * b = b * a := by
    intro a b
    obtain ⟨g1, hg1⟩ := AlgEquiv.restrictNormalHom_surjective (F := K.carrier)
      (K₁ := (Lt : Type _)) K.closure a
    obtain ⟨g2, hg2⟩ := AlgEquiv.restrictNormalHom_surjective (F := K.carrier)
      (K₁ := (Lt : Type _)) K.closure b
    have hk : g1 * g2 * g1⁻¹ * g2⁻¹ ∈ (Lt.fixingSubgroup : Subgroup K.absGal) := by
      rw [hfixLt]
      exact le_sup_right (a := Subgroup.zpowers σ) (hcommM g1 g2)
    rw [← ker_restrictNormalHom_eq_fixingSubgroup Lt, MonoidHom.mem_ker] at hk
    simp only [map_mul, map_inv] at hk
    rw [← hg1, ← hg2]
    exact mul_comm_of_commutator_eq_one hk
  -- ★強化点:`S ≤ Gal(K̄/F)` ゆえ `F ≤ L_t`
  have hFS : S ≤ (F.fixingSubgroup : Subgroup K.absGal) := by
    rw [hSdef]
    refine sup_le ?_ ?_
    · rw [Subgroup.zpowers_le]; exact hσF
    · rw [← hMfix]; exact IntermediateField.fixingSubgroup_le (le_trans hFL le_sup_left)
  have hFLt : F ≤ Lt := by
    intro x hx
    rw [hLtdef, IntermediateField.mem_fixedField_iff]
    intro g hg
    exact (IntermediateField.mem_fixingSubgroup_iff F g).mp (hFS hg) x hx
  -- 原始元定理で `L_t` を単項化する
  haveI : FiniteDimensional K.carrier (L ⊔ Km : IntermediateField K.carrier K.closure) :=
    IntermediateField.finiteDimensional_sup L Km
  haveI : FiniteDimensional K.carrier Lt :=
    Module.Finite.of_injective (IntermediateField.inclusion hLtle).toLinearMap
      (IntermediateField.inclusion hLtle).injective
  haveI : Algebra.IsSeparable K.carrier Lt := IntermediateField.isSeparable_tower_bot K.carrier Lt
  obtain ⟨α0, hα0⟩ := Field.exists_primitive_element K.carrier Lt
  have hadjα : IntermediateField.adjoin K.carrier ({(α0 : K.closure)} : Set K.closure) = Lt := by
    have h := IntermediateField.lift_adjoin K.carrier Lt ({α0} : Set Lt)
    rw [hα0, IntermediateField.lift_top] at h
    simpa using h.symm
  refine ⟨(α0 : K.closure), w, ?_, ?_, ?_, ?_, huw, hrankw, ?_⟩
  · exact isTotallyRamifiedAdjoin_of_inf_unramifiedClosure_eq_bot K _ (by rw [hadjα]; exact hLtur)
  · rw [hadjα]; infer_instance
  · rw [hadjα]; exact hLtcomm
  · rw [hadjα]; exact hFLt
  · rw [hadjα]; exact hsup

/-- **`unramLevel` 版**——不分岐側を `UnramifiedZhat.lean` の段 `unramLevel K m` で
名指しした形(`exists_totallyRamified_abelian_split_unramLevel` の強化版)。 -/
theorem exists_totallyRamified_abelian_split_ge_unramLevel (K : PAdicLocalField p)
    (L : IntermediateField K.carrier K.closure) [FiniteDimensional K.carrier L]
    [Normal K.carrier L] (hab : ∀ a b : (L ≃ₐ[K.carrier] L), a * b = b * a)
    {m : ℕ} (hm : m ≠ 0)
    (hexp : Monoid.exponent (L ≃ₐ[K.carrier] L) ∣ m)
    (F : IntermediateField K.carrier K.closure) (hFL : F ≤ L)
    (hFur : F ⊓ unramifiedClosure K = ⊥) :
    ∃ α : K.closure,
      IsTotallyRamifiedAdjoin K α ∧
      Normal K.carrier (IntermediateField.adjoin K.carrier ({α} : Set K.closure)) ∧
      (∀ a b : ((IntermediateField.adjoin K.carrier ({α} : Set K.closure))
        ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({α} : Set K.closure))),
          a * b = b * a) ∧
      F ≤ IntermediateField.adjoin K.carrier ({α} : Set K.closure) ∧
      IntermediateField.adjoin K.carrier ({α} : Set K.closure)
          ⊔ IntermediateField.lift (unramLevel K m)
        = L ⊔ IntermediateField.lift (unramLevel K m) := by
  obtain ⟨α, β, htr, hnor, hcomm, hFα, huβ, hrβ, heq⟩ :=
    exists_totallyRamified_abelian_split_ge K L hab hm hexp F hFL hFur
  have hβ : IntermediateField.adjoin K.carrier ({β} : Set K.closure)
      = IntermediateField.lift (unramLevel K m) := by
    rw [lift_unramLevel]
    exact adjoin_eq_of_isUnramified K β _ huβ (unramLevel_spec K hm).1
      (by rw [hrβ, (unramLevel_spec K hm).2.1])
  exact ⟨α, htr, hnor, hcomm, hFα, by rw [← hβ]; exact heq⟩

end ABC3.Found.PGC
