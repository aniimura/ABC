import ABC3.Found.PGC.DworkMultiplicative
import ABC3.Found.PGC.LubinTateTowerFIndependent
import ABC3.Found.PGC.LubinTateFieldFIndependent

/-!
# `K^m = K^ur·K(Λ_{f,m})` は素元 `π` の取り方に依らない —— 完備化から降りる段

経路 Λ の節点 Λ6 の締め。`Found/PGC/LubinTateTowerFIndependent.lean` が
**完備化のレベル**で `K̂^m_f = K̂^m_g`(Corollary 4.9 前半)を出したので、
本ファイルはそれを**代数の側へ降ろす**。

典拠: T. Yoshida, *Local Class Field Theory via Lubin-Tate Theory*
(Ann. Fac. Sci. Toulouse 17-2, 2008; arXiv math/0606108) Definition 4.10(物理 p.9)。
構造化済み原文は `ResearchPaper/1_Structured/Local Class Field Theory via Lubin-Tate Theory/`
の `section-4.html` の `#def-4-10`。

原文 (Yoshida08 p.9):
> Definition 4.10. For any f ∈ O[scr]_L[X] with L/K finite, set K^m := K^urL^m_f. Then K^m/K is finitely ramified, and Galois by Proposition 4.7(i). By Lemma 2.2, the completion of K^m is KL^m_f = K^m and K^m = K^m ∩ K^sep, thus independent of f. Setting K^LT := _m≥1 K^m = K^LT ∩K^sep, we have W(K^LT/K) ∼ = W(K^LT/K) by the remark after Definition 2.5. We call a finite extension of K a Lubin-Tate extension if it is contained in K^LT. We call the inverse of ρ the Artin map of K and write Art_K : K^× ∼ =−→ W(K^LT/K). We have v ◦ Art_K = v.

原文 (Yoshida08 p.9):
> Corollary 4.9. The K^m_f and ρ_f,m, hence also K^LT_f and ρ_f, of Proposition 4.7(ii) do not depend on f. (We will drop the subscript f and write K^m, ρ_m, K^LT and ρ.)

## ★★★ 何が真で、何が偽か(この節点でいちばん大事な区別)

原典 Definition 4.10 が独立性を主張しているのは

  `K^m := K^ur · L^m_f`   (不分岐閉包を掛けた体)

であって、Lubin-Tate 部分だけの `K_π := K(Λ_{f,∞})` **ではない**。

★★`K_π = K_{π′}` は **偽** である。`K = ℚ_p`(`p` 奇素数)、`π = p`、`π′ = -p`
で反例になる: `K_p = ℚ_p(µ_{p^∞})` だが `Art(-1)` は `µ_{p^∞}` 上で反転として
働くので `Art(-p)` は `ℚ_p(µ_{p^∞})` を各点固定せず、`K_{-p} ≠ K_p`。
本ファイルはしたがって

  `K_π ⊔ K^ur = K_{π′} ⊔ K^ur`

の形でだけ主張する。★この `⊔ K^ur` を落とすと主張は偽になるので、
後続ノードは落とさないこと。木の側の裏付けとしては
`AbelianDecomposition.lean::lubinTateClosure_inf_unramifiedClosure`
(`K_π ⊓ K^ur = ⊥`)がある——`K_π` 自身は `K^ur` の情報を一切持たない。

## 原典との対応

原典は Lemma 2.2(iii)(`Ê ∩ K^sep = E`、henselian 体上の Krasner の補題)で
完備側の一意性を代数側へ降ろす。★本ファイルは **Lemma 2.2(iii) を経由しない**。
代わりに無限次 Galois 対応を使う:

1. `σ ∈ Gal(K^al/K)` はスペクトルノルムに関して等長(`norm_absGal`、既に木にある)
   なので、`ℂ_K = K^al^` へ一意に連続延長する(`absGalCompletionRingEquiv`)。
2. `σ` が `K^ur` を各点固定するなら、延長は稠密性で `K̂^ur` を各点固定する
   (`absGalCompletionRingEquiv_algebraMap`)——ここで初めて `ℂ_K` の中に
   居ることを使う。これで延長は `K̂^ur`-代数同型になる。
3. `σ` がさらに `Λ_{f,n}` を各点固定するなら、延長は
   `K̂^m_f = K̂^ur(Λ_{f,n})` を各点固定する(生成元で決まるから)。
4. `x ∈ Λ_{g,n}` なら `ι(x) ∈ K̂^m_g = K̂^m_f`(Corollary 4.9 前半)なので
   `ι(σx) = ι(x)`、`ι` は単射だから `σx = x`。
5. `σ` は `K^ur ⊔ K(Λ_{f,n})` を固定する任意の元だったので、無限次 Galois 対応
   (`InfiniteGalois.fixedField_fixingSubgroup`)より `x ∈ K^ur ⊔ K(Λ_{f,n})`。

★段 2 が「どこで `K̄`(正しくは `ℂ_K`)の中に居ることを使ったか」である:
2 つの体を**同じ完備体 `ℂ_K` の中**で見ているので、`K^ur` の完備化 `K̂^ur` と
`Λ` の生成する体が 1 つの体の中で会う。これを落とすと「同じ完備化を持つ
別々の体」を排除できない。

## 退化の自己検査

* `ϖ = uπ`(`u` 単数)を落とすと偽。★任意の 2 元では駄目で、たとえば `π` と
  `π²` では `Ideal.span {π²} ≠ maximalIdeal` なので `lubinTateCompletionField`
  の入力条件(`maximalIdeal = span {ϖ}`)がそもそも作れない。逆に言えば、
  「素元どうし」という条件は statement の `hϖmax` に埋め込まれている
  (`lubinTateClosure_sup_unramifiedClosure_eq_of_uniformizer` はそこから
  `Associated π ϖ` を取り出すだけ)。
* `⊔ unramifiedClosure K` を落とすと偽(上の反例)。
* Dwork(Λ6)の `σ` は `∃ σ, ∀ u, ∃ ξ` の内側にある。本ファイルは Dwork を
  **直接は呼ばない**——`lubinTateCompletionField_eq` の中で消費済みで、
  そこから出てくるのは体の等式だけなので、量化の順は壊しようがない。

## 逸脱の記録

1. 原典 Definition 4.10 は `f ∈ 𝒪_L[[X]]`(`L/K` 有限不分岐)を許すが、本ファイルは
   **`L = K`**、すなわち `f ∈ 𝒪_K[[X]]` の場合だけを扱う。木の Lubin-Tate 機械
   (`lubinTateClosure` / `lubinTateCompletionField`)が `𝒪_K` 係数で建っているため。
   `K^ur·K(Λ_{f,m})` は `L` を変えても変わらないので、後続の主張には影響しない。
2. 原典は Lemma 2.2(iii) を経由するが、本ファイルは無限次 Galois 対応で降ろす
   (上記)。★どちらも結論は同じ等式で、後続が消費するのは等式だけ。
   ★Lemma 2.2(iii) 自身(`Ê ∩ K^sep = E`)は本ファイルでは**証明していない**。
3. 原典は `ρ_{f,m} = ρ_{f′,m}`(写像の一致)も同じ Corollary で述べるが、
   それは Lemma 4.5 を要するので本ファイルには含まれない(別ノード)。

## 抽象核と具体層

| 段 | 宣言 | 分岐・付値・Lubin-Tate の語彙 |
|---|---|---|
| 抽象核 A | `fixedIntermediateField` / `algEquiv_eq_self_of_mem_adjoin` | 出てこない(体と `AlgEquiv` だけ) |
| 抽象核 B | `mem_of_forall_fixingSubgroup_eq` | 出てこない(無限次 Galois だけ) |
| 具体層 1 | `absGalCompletionRingEquiv` ほか | `ℂ_K` |
| 具体層 2 | `absGalCompletion_eq_self_of_mem_lubinTateCompletionField` | Lubin-Tate |
| 具体層 3 | `lubinTateLevelField_le_sup_of_completionField_eq` | 降下の本体 |
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NNReal Valued Classical

variable {p : ℕ} [Fact p.Prime]

def lubinTateLevelField_sup_unramifiedClosure_eq_of_uniformizer.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 9, item := "Definition 4.10", sectionId := "def-4-10" }

def lubinTateClosure_sup_unramifiedClosure_eq_of_uniformizer.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 9, item := "Definition 4.10", sectionId := "def-4-10" }

/-! ★以下の 2 本は Corollary 4.9(完備化のレベルの独立性)を降下と繋いだ形なので、
`.src` は Corollary 4.9 を指す。 -/

def lubinTateLevelField_sup_unramifiedClosure_eq_of_isUnit_mul.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 9, item := "Corollary 4.9", sectionId := "cor-4-9" }

def lubinTateClosure_sup_unramifiedClosure_eq_of_isUnit_mul.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 9, item := "Corollary 4.9", sectionId := "cor-4-9" }

/-! ## 1. 抽象核 A —— 1 つの自己同型の固定体

★分岐・付値・完備化の語彙が 1 つも出てこない。体と `AlgEquiv` だけ。
mathlib の `IntermediateField.fixedField` は**部分群**を取るので、
「1 つの `φ`」で使うには `Subgroup.closure` の帰納法が要る。ここでは
1 元版を直に作って、`IntermediateField.adjoin_le_iff` に噛ませる。 -/

/-- `φ : Ω ≃ₐ[k] Ω` の固定体を `IntermediateField k Ω` として。 -/
def fixedIntermediateField {k Ω : Type*} [Field k] [Field Ω] [Algebra k Ω]
    (φ : Ω ≃ₐ[k] Ω) : IntermediateField k Ω where
  carrier := {x | φ x = x}
  mul_mem' := by
    intro a b ha hb
    simp only [Set.mem_setOf_eq] at ha hb ⊢
    rw [map_mul, ha, hb]
  one_mem' := by simp [Set.mem_setOf_eq]
  add_mem' := by
    intro a b ha hb
    simp only [Set.mem_setOf_eq] at ha hb ⊢
    rw [map_add, ha, hb]
  zero_mem' := by simp [Set.mem_setOf_eq]
  algebraMap_mem' := fun r => φ.commutes r
  inv_mem' := by
    intro a ha
    simp only [Set.mem_setOf_eq] at ha ⊢
    rw [map_inv₀, ha]

@[simp] theorem mem_fixedIntermediateField_iff {k Ω : Type*} [Field k] [Field Ω] [Algebra k Ω]
    (φ : Ω ≃ₐ[k] Ω) (x : Ω) : x ∈ fixedIntermediateField φ ↔ φ x = x := Iff.rfl

/-- **生成元を固定する `k`-自己同型は `adjoin k S` を各点固定する。**

★これが「`σ` が `K̂^ur` と `Λ_{f,n}` を固定するなら `K̂^m_f` を固定する」の中身で、
分岐の言葉は 1 つも要らない。 -/
theorem algEquiv_eq_self_of_mem_adjoin {k Ω : Type*} [Field k] [Field Ω] [Algebra k Ω]
    (φ : Ω ≃ₐ[k] Ω) (S : Set Ω) (hS : ∀ s ∈ S, φ s = s)
    {x : Ω} (hx : x ∈ IntermediateField.adjoin k S) : φ x = x :=
  (mem_fixedIntermediateField_iff φ x).mp
    ((IntermediateField.adjoin_le_iff.mpr
      (fun s hs => (mem_fixedIntermediateField_iff φ s).mpr (hS s hs))) hx)

/-! ## 2. 抽象核 B —— 無限次 Galois 対応の「使う向き」

`Ω/k` が(無限次でよい)Galois なら、`E` を各点固定する自己同型がすべて `x` を
固定するとき `x ∈ E`。★mathlib の `InfiniteGalois.fixedField_fixingSubgroup` を
使いやすい形に言い換えただけだが、これを名前で切り出しておくと降下の段が
1 行になる。 -/

/-- **`E` を固定する自己同型がすべて `x` を固定するなら `x ∈ E`。** -/
theorem mem_of_forall_fixingSubgroup_eq {k Ω : Type*} [Field k] [Field Ω] [Algebra k Ω]
    [IsGalois k Ω] (E : IntermediateField k Ω) (x : Ω)
    (h : ∀ σ : Ω ≃ₐ[k] Ω, (∀ y ∈ E, σ y = y) → σ x = x) : x ∈ E := by
  rw [← InfiniteGalois.fixedField_fixingSubgroup E]
  exact fun σ => h σ.1 ((IntermediateField.mem_fixingSubgroup_iff E σ.1).mp σ.2)

/-! ## 3. 具体層 1 —— `Gal(K^al/K)` の `ℂ_K` への連続延長

`UnramifiedCompletion.lean` が `Gal(K^ur/K) → Aut(K̂^{ur})` でやったことを、
`K^al` と `ℂ_K` で繰り返す。等長性は木に既にある `norm_absGal` から。 -/

/-- `σ ∈ Gal(K^al/K)` は `K^al` 上等長(`norm_absGal` の言い換え)。 -/
theorem isometry_absGal (K : PAdicLocalField p) (σ : K.absGal) :
    Isometry (σ : K.closure → K.closure) :=
  Isometry.of_dist_eq fun x y => by
    rw [dist_eq_norm, dist_eq_norm, ← map_sub, norm_absGal]

/-- `σ ∈ Gal(K^al/K)` の `ℂ_K` への延長(環同型として)。 -/
noncomputable def absGalCompletionRingEquiv (K : PAdicLocalField p) (σ : K.absGal) :
    closureCompletion K ≃+* closureCompletion K :=
  UniformSpace.Completion.mapRingEquiv σ.toRingEquiv (isometry_absGal K σ).continuous
    (isometry_absGal K σ.symm).continuous

@[simp] theorem absGalCompletionRingEquiv_coe (K : PAdicLocalField p) (σ : K.absGal)
    (x : K.closure) :
    absGalCompletionRingEquiv K σ (x : closureCompletion K)
      = ((σ x : K.closure) : closureCompletion K) :=
  UniformSpace.Completion.map_coe (isometry_absGal K σ).uniformContinuous x

theorem continuous_absGalCompletionRingEquiv (K : PAdicLocalField p) (σ : K.absGal) :
    Continuous (absGalCompletionRingEquiv K σ) := UniformSpace.Completion.continuous_map

/-- ★★**`σ` が `K^ur` を各点固定するなら、延長は `K̂^{ur}` を各点固定する。**

`K^ur` は `K̂^{ur}` で稠密(`denseRange_coe_unramifiedCompletion`)で、
両辺とも連続だから `DenseRange.equalizer` で一致する。

★ここが「完備化から降りる」ために `ℂ_K` の中に居ることを使う唯一の場所である
——`K^ur` の完備化 `K̂^{ur}` と `K^al` の元が、同じ体 `ℂ_K` の中で会っている。 -/
theorem absGalCompletionRingEquiv_algebraMap (K : PAdicLocalField p) (σ : K.absGal)
    (hσ : ∀ y ∈ unramifiedClosure K, σ y = y) (z : unramifiedCompletion K) :
    absGalCompletionRingEquiv K σ
        (algebraMap (unramifiedCompletion K) (closureCompletion K) z)
      = algebraMap (unramifiedCompletion K) (closureCompletion K) z := by
  have hfun : (fun w : unramifiedCompletion K => absGalCompletionRingEquiv K σ
        (algebraMap (unramifiedCompletion K) (closureCompletion K) w))
      = (fun w : unramifiedCompletion K =>
        algebraMap (unramifiedCompletion K) (closureCompletion K) w) := by
    refine (denseRange_coe_unramifiedCompletion K).equalizer
      ((continuous_absGalCompletionRingEquiv K σ).comp
        (continuous_unramifiedToClosureCompletion K))
      (continuous_unramifiedToClosureCompletion K) ?_
    funext y
    show absGalCompletionRingEquiv K σ
        (unramifiedToClosureCompletion K (y : unramifiedCompletion K))
      = unramifiedToClosureCompletion K (y : unramifiedCompletion K)
    rw [unramifiedToClosureCompletion_coe, absGalCompletionRingEquiv_coe, hσ (y : K.closure) y.2]
  exact congrFun hfun z

/-- `K^ur` を各点固定する `σ` の延長を **`K̂^{ur}`-代数同型**として見る。 -/
noncomputable def absGalCompletionAlg (K : PAdicLocalField p) (σ : K.absGal)
    (hσ : ∀ y ∈ unramifiedClosure K, σ y = y) :
    closureCompletion K ≃ₐ[unramifiedCompletion K] closureCompletion K :=
  AlgEquiv.ofRingEquiv (f := absGalCompletionRingEquiv K σ)
    (absGalCompletionRingEquiv_algebraMap K σ hσ)

@[simp] theorem absGalCompletionAlg_apply (K : PAdicLocalField p) (σ : K.absGal)
    (hσ : ∀ y ∈ unramifiedClosure K, σ y = y) (z : closureCompletion K) :
    absGalCompletionAlg K σ hσ z = absGalCompletionRingEquiv K σ z := rfl

/-! ## 4. 具体層 2 —— 延長は `K̂^m_f` を各点固定する

`K̂^m_f = K̂^{ur}(ι(Λ_{f,n}))` は定義から `adjoin` なので、抽象核 A に
「生成元 = `Λ_{f,n}` の像」を渡すだけ。 -/

/-- **`σ` が `K^ur` と `Λ_{f,n}` を各点固定するなら、延長は `K̂^m_f` を各点固定する。** -/
theorem absGalCompletion_eq_self_of_mem_lubinTateCompletionField
    (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (σ : K.absGal)
    (hσur : ∀ y ∈ unramifiedClosure K, σ y = y)
    (hσΛ : ∀ y ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n, σ y = y)
    {z : closureCompletion K}
    (hz : z ∈ lubinTateCompletionField K hq hπmax hπne0 f hf0 hf1 hf n) :
    absGalCompletionRingEquiv K σ z = z := by
  rw [← absGalCompletionAlg_apply K σ hσur]
  refine algEquiv_eq_self_of_mem_adjoin (absGalCompletionAlg K σ hσur) _ ?_ hz
  rintro s ⟨lam, hlam, rfl⟩
  rw [absGalCompletionAlg_apply, closureCompletionCoe_apply, absGalCompletionRingEquiv_coe,
    hσΛ lam (Finset.mem_coe.mp hlam)]

/-! ## 5. `K^m` の分岐部分 `K(Λ_{f,m})` -/

/-- **`L^m_f` の代数版** —— `K(Λ_{f,n}) ⊆ K^al`。原典 Definition 4.10 の
`K^m := K^ur L^m_f` の第 2 因子(ただし原典の `L` を `K` に取った場合、逸脱 1)。 -/
noncomputable def lubinTateLevelField (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) : IntermediateField K.carrier K.closure :=
  IntermediateField.adjoin K.carrier
    ((iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n : Finset K.closure) :
      Set K.closure)

/-- `K(Λ_{f,n}) ≤ K_π`。 -/
theorem lubinTateLevelField_le_lubinTateClosure (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) :
    lubinTateLevelField K hq hπmax hπne0 f hf0 hf1 hf n
      ≤ lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf := by
  rw [lubinTateLevelField, IntermediateField.adjoin_le_iff]
  intro x hx
  refine IntermediateField.subset_adjoin K.carrier _ ?_
  simp only [lubinTateTorsionSet, Set.mem_iUnion]
  exact ⟨n, hx⟩

/-! ## 6. 具体層 3 —— 降下(completion から代数へ)

★本ファイルの本体。`K̂^m_g = K̂^m_f` から `K(Λ_{g,n}) ≤ K(Λ_{f,n})·K^ur` を出す。 -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**降下の段(段数固定)** —— `K̂^m_g = K̂^m_f` なら
`K(Λ_{g,n}) ≤ K(Λ_{f,n}) ⊔ K^ur`。

★仮定 `hcomp` は「完備化のレベルの等式」だけ。`π` と `ϖ` の関係
(`ϖ = uπ`)は `hcomp` の**中に**入っている(それを供給するのが
`lubinTateCompletionField_eq` = Corollary 4.9 前半 = Dwork の消費側)。
したがって本補題は左右対称に使える。

証明は無限次 Galois 対応 —— `E := K(Λ_{f,n}) ⊔ K^ur` を各点固定する
`σ ∈ Gal(K^al/K)` を取り、`ℂ_K` へ延長して `K̂^m_f = K̂^m_g ∋ ι(x)` を
固定させ、`ι` の単射性で `σ x = x` に戻す。 -/
theorem lubinTateLevelField_le_sup_of_completionField_eq
    (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    {ϖ : 𝒪[K.carrier]} (hϖmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {ϖ})
    (hϖne0 : ϖ ≠ 0)
    (g : PowerSeries (𝒪[K.carrier])) (hg0 : PowerSeries.coeff 0 g = 0)
    (hg1 : PowerSeries.coeff 1 g = ϖ)
    (hg : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) g = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ)
    (hcomp : lubinTateCompletionField K hq hϖmax hϖne0 g hg0 hg1 hg n
        = lubinTateCompletionField K hq hπmax hπne0 f hf0 hf1 hf n) :
    lubinTateLevelField K hq hϖmax hϖne0 g hg0 hg1 hg n
      ≤ lubinTateLevelField K hq hπmax hπne0 f hf0 hf1 hf n ⊔ unramifiedClosure K := by
  haveI := isGalois_closure K
  rw [lubinTateLevelField, IntermediateField.adjoin_le_iff]
  intro x hx
  refine mem_of_forall_fixingSubgroup_eq _ x ?_
  intro σ hσ
  have hσur : ∀ y ∈ unramifiedClosure K, σ y = y :=
    fun y hy => hσ y (SetLike.le_def.mp le_sup_right hy)
  have hσΛ : ∀ y ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n, σ y = y := by
    intro y hy
    exact hσ y (SetLike.le_def.mp le_sup_left
      (IntermediateField.subset_adjoin K.carrier _ (Finset.mem_coe.mpr hy)))
  have hmem : ((x : K.closure) : closureCompletion K)
      ∈ lubinTateCompletionField K hq hπmax hπne0 f hf0 hf1 hf n := by
    rw [← hcomp]
    exact IntermediateField.subset_adjoin _ _ ⟨x, hx, rfl⟩
  have hfix := absGalCompletion_eq_self_of_mem_lubinTateCompletionField K hq hπmax hπne0
    f hf0 hf1 hf n σ hσur hσΛ hmem
  rw [absGalCompletionRingEquiv_coe] at hfix
  exact closureCompletionCoe_injective K hfix

/-- **降下の段(合併)** —— 全段で `K̂^m_g = K̂^m_f` なら `K_{ϖ,g} ≤ K_{π,f} ⊔ K^ur`。 -/
theorem lubinTateClosure_le_sup_of_completionField_eq
    (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    {ϖ : 𝒪[K.carrier]} (hϖmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {ϖ})
    (hϖne0 : ϖ ≠ 0)
    (g : PowerSeries (𝒪[K.carrier])) (hg0 : PowerSeries.coeff 0 g = 0)
    (hg1 : PowerSeries.coeff 1 g = ϖ)
    (hg : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) g = PowerSeries.X ^ (pp ^ ff))
    (hcomp : ∀ n : ℕ, lubinTateCompletionField K hq hϖmax hϖne0 g hg0 hg1 hg n
        = lubinTateCompletionField K hq hπmax hπne0 f hf0 hf1 hf n) :
    lubinTateClosure K hq hϖmax hϖne0 g hg0 hg1 hg
      ≤ lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ⊔ unramifiedClosure K := by
  rw [lubinTateClosure, IntermediateField.adjoin_le_iff]
  intro x hx
  simp only [lubinTateTorsionSet, Set.mem_iUnion, Finset.mem_coe] at hx
  obtain ⟨n, hxn⟩ := hx
  have hlevel := lubinTateLevelField_le_sup_of_completionField_eq K hq hπmax hπne0 f hf0 hf1 hf
    hϖmax hϖne0 g hg0 hg1 hg n (hcomp n)
  have hx' : x ∈ lubinTateLevelField K hq hϖmax hϖne0 g hg0 hg1 hg n :=
    IntermediateField.subset_adjoin K.carrier _ (Finset.mem_coe.mpr hxn)
  refine SetLike.le_def.mp (sup_le_sup_right ?_ _) (SetLike.le_def.mp hlevel hx')
  exact lubinTateLevelField_le_lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf n

/-! ## 7. 主定理 —— `K^m` と `K^LT` は素元の取り方に依らない -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`K^m` は素元に依らない(単数倍の形)** ——
`f ∈ F_π`、`g ∈ F_{uπ}`(`u` 単数)に対し

`K^ur · K(Λ_{f,n}) = K^ur · K(Λ_{g,n})`。

これが原典 Definition 4.10 の「thus independent of f」の体の部分である。

★★`⊔ unramifiedClosure K` を落とすと**偽**(モジュール docstring の反例)。 -/
theorem lubinTateLevelField_sup_unramifiedClosure_eq_of_isUnit_mul
    (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    {u : 𝒪[K.carrier]} (hu : IsUnit u)
    (hϖmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {u * π})
    (hϖne0 : u * π ≠ 0)
    (g : PowerSeries (𝒪[K.carrier])) (hg0 : PowerSeries.coeff 0 g = 0)
    (hg1 : PowerSeries.coeff 1 g = u * π)
    (hg : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) g = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) :
    lubinTateLevelField K hq hπmax hπne0 f hf0 hf1 hf n ⊔ unramifiedClosure K
      = lubinTateLevelField K hq hϖmax hϖne0 g hg0 hg1 hg n ⊔ unramifiedClosure K := by
  have hcomp : lubinTateCompletionField K hq hϖmax hϖne0 g hg0 hg1 hg n
      = lubinTateCompletionField K hq hπmax hπne0 f hf0 hf1 hf n :=
    (lubinTateCompletionField_eq K hq hπmax hπne0 f hf0 hf1 hf hu g hg0 hg1 hg n).symm
  refine le_antisymm (sup_le ?_ le_sup_right) (sup_le ?_ le_sup_right)
  · exact lubinTateLevelField_le_sup_of_completionField_eq K hq hϖmax hϖne0 g hg0 hg1 hg
      hπmax hπne0 f hf0 hf1 hf n hcomp.symm
  · exact lubinTateLevelField_le_sup_of_completionField_eq K hq hπmax hπne0 f hf0 hf1 hf
      hϖmax hϖne0 g hg0 hg1 hg n hcomp

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`K^LT` は素元に依らない(単数倍の形)** ——
`K^ur · K_{π,f} = K^ur · K_{uπ,g}`。 -/
theorem lubinTateClosure_sup_unramifiedClosure_eq_of_isUnit_mul
    (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    {u : 𝒪[K.carrier]} (hu : IsUnit u)
    (hϖmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {u * π})
    (hϖne0 : u * π ≠ 0)
    (g : PowerSeries (𝒪[K.carrier])) (hg0 : PowerSeries.coeff 0 g = 0)
    (hg1 : PowerSeries.coeff 1 g = u * π)
    (hg : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) g = PowerSeries.X ^ (pp ^ ff)) :
    lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ⊔ unramifiedClosure K
      = lubinTateClosure K hq hϖmax hϖne0 g hg0 hg1 hg ⊔ unramifiedClosure K := by
  have hcomp : ∀ n : ℕ, lubinTateCompletionField K hq hϖmax hϖne0 g hg0 hg1 hg n
      = lubinTateCompletionField K hq hπmax hπne0 f hf0 hf1 hf n := fun n =>
    (lubinTateCompletionField_eq K hq hπmax hπne0 f hf0 hf1 hf hu g hg0 hg1 hg n).symm
  refine le_antisymm (sup_le ?_ le_sup_right) (sup_le ?_ le_sup_right)
  · exact lubinTateClosure_le_sup_of_completionField_eq K hq hϖmax hϖne0 g hg0 hg1 hg
      hπmax hπne0 f hf0 hf1 hf (fun n => (hcomp n).symm)
  · exact lubinTateClosure_le_sup_of_completionField_eq K hq hπmax hπne0 f hf0 hf1 hf
      hϖmax hϖne0 g hg0 hg1 hg hcomp

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`K^m` は素元に依らない(任意の 2 素元)** ——
`π`・`ϖ` を `K` の任意の素元、`f ∈ F_π`・`g ∈ F_ϖ` とすると

`K^ur · K(Λ_{f,n}) = K^ur · K(Λ_{g,n})`。

★2 つの素元は自動的に単数倍で結ばれる(`Ideal.span_singleton_eq_span_singleton`
が `Associated π ϖ` を返す)ので、単数倍の形からそのまま出る。 -/
theorem lubinTateLevelField_sup_unramifiedClosure_eq_of_uniformizer
    (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    {ϖ : 𝒪[K.carrier]} (hϖmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {ϖ})
    (hϖne0 : ϖ ≠ 0)
    (g : PowerSeries (𝒪[K.carrier])) (hg0 : PowerSeries.coeff 0 g = 0)
    (hg1 : PowerSeries.coeff 1 g = ϖ)
    (hg : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) g = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) :
    lubinTateLevelField K hq hπmax hπne0 f hf0 hf1 hf n ⊔ unramifiedClosure K
      = lubinTateLevelField K hq hϖmax hϖne0 g hg0 hg1 hg n ⊔ unramifiedClosure K := by
  obtain ⟨u, hu⟩ : Associated π ϖ :=
    Ideal.span_singleton_eq_span_singleton.mp (hπmax.symm.trans hϖmax)
  have hϖ : ϖ = (u : 𝒪[K.carrier]) * π := by rw [← hu]; ring
  subst hϖ
  exact lubinTateLevelField_sup_unramifiedClosure_eq_of_isUnit_mul K hq hπmax hπne0 f hf0 hf1 hf
    u.isUnit hϖmax hϖne0 g hg0 hg1 hg n

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`K^LT` は素元に依らない(任意の 2 素元)** ——
`π`・`ϖ` を `K` の任意の素元、`f ∈ F_π`・`g ∈ F_ϖ` とすると

`K^ur · K_{π,f} = K^ur · K_{ϖ,g}`。

★★これが原典 Corollary 4.9 + Definition 4.10 の到達点の代数側である。
★★`K_{π,f} = K_{ϖ,g}` は偽なので、`⊔ unramifiedClosure K` は落とせない。 -/
theorem lubinTateClosure_sup_unramifiedClosure_eq_of_uniformizer
    (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    {ϖ : 𝒪[K.carrier]} (hϖmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {ϖ})
    (hϖne0 : ϖ ≠ 0)
    (g : PowerSeries (𝒪[K.carrier])) (hg0 : PowerSeries.coeff 0 g = 0)
    (hg1 : PowerSeries.coeff 1 g = ϖ)
    (hg : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) g = PowerSeries.X ^ (pp ^ ff)) :
    lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ⊔ unramifiedClosure K
      = lubinTateClosure K hq hϖmax hϖne0 g hg0 hg1 hg ⊔ unramifiedClosure K := by
  obtain ⟨u, hu⟩ : Associated π ϖ :=
    Ideal.span_singleton_eq_span_singleton.mp (hπmax.symm.trans hϖmax)
  have hϖ : ϖ = (u : 𝒪[K.carrier]) * π := by rw [← hu]; ring
  subst hϖ
  exact lubinTateClosure_sup_unramifiedClosure_eq_of_isUnit_mul K hq hπmax hπne0 f hf0 hf1 hf
    u.isUnit hϖmax hϖne0 g hg0 hg1 hg

end ABC3.Found.PGC
