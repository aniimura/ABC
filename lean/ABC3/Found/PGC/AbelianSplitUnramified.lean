import ABC3.Found.PGC.AbelianFrobeniusSplit
import ABC3.Found.PGC.UnramifiedZhat

/-!
# Milne Lemma 4.11 と Lemma 4.10——有限アーベル拡大を「全分岐 × 不分岐」に分ける

Milne, *Class Field Theory*, p.57 の Lemma 4.10 と Lemma 4.11。
経路 Λ7(`K^ab = K_π · K^ur`)の起点である。

## 原典

Lemma 4.10(p.57):K の有限不分岐拡大はすべて `K_π · K^ur` に含まれる。

Lemma 4.11(p.57):L を K の指数 m の有限アーベル拡大、K_m を K の次数 m の
不分岐拡大とするとき、K 上完全分岐なアーベル拡大 L_t が在って
`L · K_m = L_t · K_m` となる。

## 主定理

* `exists_totallyRamified_abelian_split`——Lemma 4.11。
  `L/K` が有限アーベルで `Monoid.exponent Gal(L/K) ∣ m`(`m ≠ 0`)のとき、
  `α β : K.closure` が在って
  - `K(α)/K` は完全分岐で、しかも `K` 上正規かつ Galois 群が可換
  - `K(β)/K` は不分岐で次数 `m`(= `K_m`)
  - `K(α) ⊔ K(β) = L ⊔ K(β)`(これが原典の `L_t · K_m = L · K_m`)
  ★`α`・`β` はすべて `∃` の内側に閉じ込めてある。
* `exists_totallyRamified_abelian_split_unramLevel`——同じものを
  `UnramifiedZhat.lean` の `unramLevel K m` を `K.closure` へ持ち上げた形で
  述べ直したもの(`IntermediateField.lift`)。
* `adjoin_le_sup_unramifiedClosure`——Lemma 4.10。

## 途中で足した在庫(本ファイルの副産物)

* `exists_unramified_subextension`——**最大不分岐部分拡大の存在**。
  任意の `x : K.closure` に対し `K(z) ⊆ K(x)` で `K(z)/K` が不分岐、
  `[K(z):K] = f(K(x)/K)`(慣性次数)なるものが在る。
  `UnramifiedExtension.lean::exists_integral_generator` の証明の
  **前半だけ**(不分岐性の仮定を使う手前まで)を取り出したもの。
* `isTotallyRamifiedAdjoin_of_inf_unramifiedClosure_eq_bot`——
  `TotallyRamified.lean::totallyRamifiedAdjoin_inf_unramifiedClosure`
  (完全分岐 ⟹ `K(α) ⊓ K^ur = ⊥`)の**逆**。上の最大不分岐部分拡大で閉じる。

## ★退化の自己検査(原典が反例を書いてくれている珍しい例)

Milne Example 4.13(p.57–58、`ℚ_5` 上の具体例)が次の 2 つを実証している。

1. **`Monoid.exponent Gal(L/K) ∣ m` を落とすと主張は偽**。
   `m` を指数より小さく取ると `σ` の `Gal(LK_m/K)` での位数が `m` を超え、
   `⟨σ⟩ ⊓ Gal(LK_m/K_m) ≠ 1` になって固定体の交わりが `K` にならない。
   本ファイルではこの仮定は `zpowers_sup_inf_eq` の仮定
   「`σ^k ∈ H_m ⟹ σ^k ∈ H_M`」に化けており、そこが崩れる。
2. **`K_m` で拡大せず `L` のまま分解しようとすると偽**。
   Example 4.13 の `L`(`ℚ_5(ζ_5, 5^{1/4})` の中の巡回 4 次拡大)は
   `L = L_t · L_u`(`L_t/K` 完全分岐・`L_u/K` 不分岐)の形に書けない。
   したがって結論を `L = K(α) ⊔ K(β)` に強めることはできず、
   `L ⊔ K(β) = K(α) ⊔ K(β)` が最良である。

## 逸脱(記録)

* 原典は `L_t := L^{⟨σ⟩}`(`LK_m` の中の `⟨σ⟩` の固定体)と置くが、
  ここでは `Γ_K` の中の開部分群 `S := ⟨σ⟩ ⊔ Gal(K̄/LK_m)` の固定体として
  同じものを作っている(商群を作らないという設計上の制約に合わせた)。
* 原典は `Gal(LK_m/K) = ⟨σ⟩ × H` という**直積分解**を経由するが、本証明は
  分解も位数の積も使わない(`AbelianFrobeniusSplit.lean` の冒頭を参照)。
* `IsTotallyRamifiedAdjoin` は単項生成でしか定義されていないので、
  `L_t` を原始元定理(`Field.exists_primitive_element`)で単項化してから渡している。
  ここで `FiniteDimensional K.carrier L_t` が要り、`L_t ≤ L ⊔ K_m` と
  `IntermediateField.finiteDimensional_sup` から取っている。
* **アーベル性の仮定 `hab` は結論のうち「`Gal(K(α)/K)` が可換」と
  「`K(α)/K` が正規」にしか使っていない**。完全分岐性と
  `K(α) ⊔ K(β) = L ⊔ K(β)` だけなら `L/K` は有限 Galois で指数が `m` を
  割りさえすればよい(原典より広い)。下流(Lemma 4.9)がアーベル性を
  要求するので、結論に含める形で残してある。
* Lemma 4.10 は本木では 1 行で済む。`unramifiedClosure K` を
  「不分岐単項拡大すべての上限」として定義しており(`UnramifiedExtension.lean`)、
  原典が `ANT 7.58` に投げている「不分岐拡大の合成はまた不分岐」の部分は
  `exists_isUnramified_ge`(不分岐拡大の有向性)が既に吸収しているため。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued

variable {p : ℕ} [Fact p.Prime]

/-! ## 0. Milne Lemma 4.10 -/

/-- **Milne Lemma 4.10**——`K` の有限不分岐拡大は `E · K^ur` に含まれる
(特に `E = K_π` で原典の形)。

本木では `unramifiedClosure K` が不分岐単項拡大すべての上限として
定義されているので 1 行。原典が `ANT 7.58` に投げている
「不分岐拡大どうしは合成できる」は `exists_isUnramified_ge` が担っている。 -/
theorem adjoin_le_sup_unramifiedClosure (K : PAdicLocalField p)
    (E : IntermediateField K.carrier K.closure) {x : K.closure} (hx : IsUnramifiedAdjoin K x) :
    IntermediateField.adjoin K.carrier ({x} : Set K.closure) ≤ E ⊔ unramifiedClosure K :=
  le_trans (adjoin_le_unramifiedClosure K hx) le_sup_right

/-! ## 1. 最大不分岐部分拡大 -/

/-- ★★★★★★**最大不分岐部分拡大の存在**——任意の `x : K.closure` に対し、
`K(z) ⊆ K(x)` で `K(z)/K` が不分岐、`[K(z):K] = f(K(x)/K)` なる `z` が在る。

剰余体 `𝓀_{K(x)}` の原始元 `t̄` の最小多項式を `𝒪_K` へ持ち上げて `f` を作り、
Hensel(`exists_root_lift_of_residue_root`)で `f` の根 `θ ∈ K(x)` を得る。
`isUnramifiedAdjoin_of_lift` により `K(θ)/K` は不分岐で
`[K(θ):K] = deg f = [𝓀_{K(x)} : 𝓀_K] = f(K(x)/K)`。

★`UnramifiedExtension.lean::exists_integral_generator` はこの構成を
「`K(x)/K` が不分岐」の仮定つきで行い `K(θ) = K(x)` まで言うが、
**前半は仮定を使っていない**。ここではその前半だけを取り出している。 -/
theorem exists_unramified_subextension (K : PAdicLocalField p) (x : K.closure) :
    ∃ z : K.closure, IsUnramifiedAdjoin K z ∧
      Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({z} : Set K.closure))
        = inertiaDegree K x ∧
      IntermediateField.adjoin K.carrier ({z} : Set K.closure)
        ≤ IntermediateField.adjoin K.carrier ({x} : Set K.closure) := by
  haveI := module_finite_adjoinIntegers K x
  haveI : Finite 𝓀[K.carrier] := residueField_finite K
  haveI : IsFractionRing 𝒪[K.carrier] K.carrier :=
    ValuationRing.instIsFractionRingInteger (K := K.carrier) Valued.v
  obtain ⟨tb, htb⟩ := Field.exists_primitive_element 𝓀[K.carrier]
    (IsLocalRing.ResidueField (adjoinIntegers K x))
  have htbint : IsIntegral 𝓀[K.carrier] tb := IsIntegral.of_finite _ _
  have hgbm : (minpoly 𝓀[K.carrier] tb).Monic := minpoly.monic htbint
  have hgbi : Irreducible (minpoly 𝓀[K.carrier] tb) := minpoly.irreducible htbint
  have hgbdeg : (minpoly 𝓀[K.carrier] tb).natDegree
      = Module.finrank 𝓀[K.carrier] (IsLocalRing.ResidueField (adjoinIntegers K x)) := by
    rw [← IntermediateField.adjoin.finrank htbint, htb, IntermediateField.finrank_top']
  have hlifts : (minpoly 𝓀[K.carrier] tb) ∈ Polynomial.lifts (IsLocalRing.residue 𝒪[K.carrier]) :=
    Polynomial.lifts_iff_coeff_lifts _ |>.mpr (fun n => ⟨_, Quotient.out_eq' _⟩)
  obtain ⟨f, hfmap, hfdeg, hfm⟩ := Polynomial.lifts_and_degree_eq_and_monic hlifts hgbm
  have hgi : Irreducible (f.map (IsLocalRing.residue 𝒪[K.carrier])) := by rw [hfmap]; exact hgbi
  have hbar : Polynomial.aeval tb (f.map (IsLocalRing.residue 𝒪[K.carrier])) = 0 := by
    rw [hfmap]; exact minpoly.aeval _ _
  obtain ⟨θA, haev⟩ := exists_root_lift_of_residue_root K x f hfm hgi tb hbar
  obtain ⟨hrankθ, huθ, -, -⟩ := isUnramifiedAdjoin_of_lift K
    (((θA : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure)) f hfm hgi haev
  have hnd : (minpoly 𝓀[K.carrier] tb).natDegree = f.natDegree := by
    rw [← hfmap, hfm.natDegree_map]
  refine ⟨((θA : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure), huθ,
    ?_, ?_⟩
  · rw [hrankθ, ← hnd, hgbdeg, ← inertiaDegree_eq_finrank_residueField K x]
  · rw [IntermediateField.adjoin_simple_le_iff]
    exact (θA : IntermediateField.adjoin K.carrier ({x} : Set K.closure)).2

/-- ★★★★★★★★**`K(α) ⊓ K^ur = K` なら `K(α)/K` は完全分岐**——
`TotallyRamified.lean::totallyRamifiedAdjoin_inf_unramifiedClosure` の逆。

最大不分岐部分拡大 `K(z) ⊆ K(α)` は `K^ur` に入るので `K(z) ⊆ K(α) ⊓ K^ur = K`、
したがって `f(K(α)/K) = [K(z):K] = 1`。 -/
theorem isTotallyRamifiedAdjoin_of_inf_unramifiedClosure_eq_bot (K : PAdicLocalField p)
    (α : K.closure)
    (h : IntermediateField.adjoin K.carrier ({α} : Set K.closure) ⊓ unramifiedClosure K = ⊥) :
    IsTotallyRamifiedAdjoin K α := by
  obtain ⟨z, huz, hrank, hle⟩ := exists_unramified_subextension K α
  have hzmem : z ∈ IntermediateField.adjoin K.carrier ({α} : Set K.closure)
      ⊓ unramifiedClosure K :=
    ⟨hle (IntermediateField.mem_adjoin_simple_self _ _),
      (mem_unramifiedClosure_iff_isUnramified K z).mpr huz⟩
  rw [h] at hzmem
  have hbot : IntermediateField.adjoin K.carrier ({z} : Set K.closure) = ⊥ :=
    le_antisymm (by rw [IntermediateField.adjoin_simple_le_iff]; exact hzmem) bot_le
  show inertiaDegree K α = 1
  rw [← hrank, hbot]
  exact IntermediateField.finrank_bot

/-! ## 2. Milne Lemma 4.11 -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**Milne Lemma 4.11**——`L/K` が有限アーベルで `Monoid.exponent Gal(L/K) ∣ m`
(`m ≠ 0`)のとき、`K.closure` の中に

* `α`:`K(α)/K` は**完全分岐**、`K` 上正規、`Gal(K(α)/K)` は可換
* `β`:`K(β)/K` は**不分岐で次数 `m`**(すなわち `K_m`)

が在って `K(α) ⊔ K(β) = L ⊔ K(β)`(原典の `L_t · K_m = L · K_m`)。

## 証明の骨格(4 段)

1. `σ` を「`Gal(K_m/K)` の生成元の `Γ_K` への持ち上げ」として取る
   (`exists_unramified_frobenius_lift`)。
2. `H_M := (L ⊔ K_m).fixingSubgroup = H_L ⊓ H_m` は開かつ正規で、
   `∀ g, g^m ∈ H_M`(`L` 側は指数の仮定、`K_m` 側は位数 `m`)。
3. `S := ⟨σ⟩ ⊔ H_M` と置くと `S ⊓ H_m = H_M`(`zpowers_sup_inf_eq`)、
   `S ⊔ H_m = ⊤`(`zpowers_sup_eq_top`)。`L_t := fixedField S` に
   無限次 Galois 対応を当てて `L_t ⊔ K_m = L ⊔ K_m`、`L_t ⊓ K_m = ⊥`。
4. `L ⊔ K_m` の中の不分岐部分は次数が `m` を割る
   (`finrank_dvd_of_isUnramified_of_pow_mem`)ので `K_m` に入る。
   よって `L_t ⊓ K^ur = ⊥`、
   `isTotallyRamifiedAdjoin_of_inf_unramifiedClosure_eq_bot` で完全分岐。

★★`m ≠ 0` と `Monoid.exponent ∣ m` はどちらも落とせない(Example 4.13)。
★★結論を `L = K(α) ⊔ K(β)` に強めることもできない(同じく Example 4.13)。 -/
theorem exists_totallyRamified_abelian_split (K : PAdicLocalField p)
    (L : IntermediateField K.carrier K.closure) [FiniteDimensional K.carrier L]
    [Normal K.carrier L] (hab : ∀ a b : (L ≃ₐ[K.carrier] L), a * b = b * a)
    {m : ℕ} (hm : m ≠ 0)
    (hexp : Monoid.exponent (L ≃ₐ[K.carrier] L) ∣ m) :
    ∃ α β : K.closure,
      IsTotallyRamifiedAdjoin K α ∧
      Normal K.carrier (IntermediateField.adjoin K.carrier ({α} : Set K.closure)) ∧
      (∀ a b : ((IntermediateField.adjoin K.carrier ({α} : Set K.closure))
        ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({α} : Set K.closure))),
          a * b = b * a) ∧
      IsUnramifiedAdjoin K β ∧
      Module.finrank K.carrier
        (IntermediateField.adjoin K.carrier ({β} : Set K.closure)) = m ∧
      IntermediateField.adjoin K.carrier ({α} : Set K.closure)
          ⊔ IntermediateField.adjoin K.carrier ({β} : Set K.closure)
        = L ⊔ IntermediateField.adjoin K.carrier ({β} : Set K.closure) := by
  haveI := isGalois_closure K
  haveI := IsAlgClosure.normal K.carrier K.closure
  obtain ⟨w, σ, huw, hrankw, hkm, htopσ, hpowm, hcommm⟩ := exists_unramified_frobenius_lift K hm
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
  -- `L_t/K` のアーベル性(★ここだけが `hab` を使う)
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
  refine ⟨(α0 : K.closure), w, ?_, ?_, ?_, huw, hrankw, ?_⟩
  · exact isTotallyRamifiedAdjoin_of_inf_unramifiedClosure_eq_bot K _ (by rw [hadjα]; exact hLtur)
  · rw [hadjα]; infer_instance
  · rw [hadjα]; exact hLtcomm
  · rw [hadjα]; exact hsup

/-! ## 3. `unramLevel` を使った言い換え

`UnramifiedZhat.lean` の `unramLevel K m` は
`IntermediateField K.carrier ↥(unramifiedClosure K)` であって
`IntermediateField K.carrier K.closure` ではない。`IntermediateField.lift`
(= `map (unramifiedClosure K).val`)で持ち上げると、上の `K(β)` と一致する。 -/

/-- `unramLevel K m` を `K.closure` へ持ち上げると `K(unramLevelGen K m)`。 -/
theorem lift_unramLevel (K : PAdicLocalField p) (m : ℕ) :
    IntermediateField.lift (unramLevel K m)
      = IntermediateField.adjoin K.carrier
        ({((unramLevelGen K m : ↥(unramifiedClosure K)) : K.closure)} : Set K.closure) :=
  map_adjoin_singleton_val K (unramLevelGen K m)

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**Milne Lemma 4.11(`unramLevel` 版)**——不分岐側を `UnramifiedZhat.lean` の
段 `unramLevel K m` で名指しした形。

`K(β)` と `lift (unramLevel K m)` はどちらも次数 `m` の不分岐拡大なので、
一意性(`adjoin_eq_of_isUnramified`)で一致する。 -/
theorem exists_totallyRamified_abelian_split_unramLevel (K : PAdicLocalField p)
    (L : IntermediateField K.carrier K.closure) [FiniteDimensional K.carrier L]
    [Normal K.carrier L] (hab : ∀ a b : (L ≃ₐ[K.carrier] L), a * b = b * a)
    {m : ℕ} (hm : m ≠ 0)
    (hexp : Monoid.exponent (L ≃ₐ[K.carrier] L) ∣ m) :
    ∃ α : K.closure,
      IsTotallyRamifiedAdjoin K α ∧
      Normal K.carrier (IntermediateField.adjoin K.carrier ({α} : Set K.closure)) ∧
      (∀ a b : ((IntermediateField.adjoin K.carrier ({α} : Set K.closure))
        ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({α} : Set K.closure))),
          a * b = b * a) ∧
      IntermediateField.adjoin K.carrier ({α} : Set K.closure)
          ⊔ IntermediateField.lift (unramLevel K m)
        = L ⊔ IntermediateField.lift (unramLevel K m) := by
  obtain ⟨α, β, htr, hnor, hcomm, huβ, hrβ, heq⟩ :=
    exists_totallyRamified_abelian_split K L hab hm hexp
  have hβ : IntermediateField.adjoin K.carrier ({β} : Set K.closure)
      = IntermediateField.lift (unramLevel K m) := by
    rw [lift_unramLevel]
    exact adjoin_eq_of_isUnramified K β _ huβ (unramLevel_spec K hm).1
      (by rw [hrβ, (unramLevel_spec K hm).2.1])
  exact ⟨α, htr, hnor, hcomm, by rw [← hβ]; exact heq⟩

/-- **`L ≤ K(α) ⊔ K_m`**(Λ7 が直接使う形)。 -/
theorem exists_totallyRamified_abelian_le_sup (K : PAdicLocalField p)
    (L : IntermediateField K.carrier K.closure) [FiniteDimensional K.carrier L]
    [Normal K.carrier L] (hab : ∀ a b : (L ≃ₐ[K.carrier] L), a * b = b * a)
    {m : ℕ} (hm : m ≠ 0)
    (hexp : Monoid.exponent (L ≃ₐ[K.carrier] L) ∣ m) :
    ∃ α : K.closure,
      IsTotallyRamifiedAdjoin K α ∧
      Normal K.carrier (IntermediateField.adjoin K.carrier ({α} : Set K.closure)) ∧
      (∀ a b : ((IntermediateField.adjoin K.carrier ({α} : Set K.closure))
        ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({α} : Set K.closure))),
          a * b = b * a) ∧
      L ≤ IntermediateField.adjoin K.carrier ({α} : Set K.closure)
        ⊔ IntermediateField.lift (unramLevel K m) := by
  obtain ⟨α, htr, hnor, hcomm, heq⟩ :=
    exists_totallyRamified_abelian_split_unramLevel K L hab hm hexp
  exact ⟨α, htr, hnor, hcomm, by rw [heq]; exact le_sup_left⟩

end ABC3.Found.PGC
