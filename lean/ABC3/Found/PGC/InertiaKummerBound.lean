import ABC3.Found.PGC.InertiaKummer
import ABC3.Found.PGC.UnramifiedGalCharCount

/-!
# (C-q) の上界 `N_n(Γ_K) ≤ n · gcd(n, q−1)` —— 仮定を外す

[pGC] Proposition 1.2 への経路 C(`ResearchPaper/pgc-goal.md`)のノード F4。

`Found/PGC/InertiaKummer.lean`(第 1033)は上界を

  `contHomCard_absGal_le_of_card_ker` ——
  `#ker(restrictInertia) ≤ n` を**仮定に残した形**

で用意していた。その一点は `Found/PGC/UnramifiedGalCharCount.lean`(第 1035、F3)の

  `contHomCard_unramifiedClosureGal : N_n(Gal(K^ur/K)) = n`

に帰着する。本ファイルはその帰着を実行し、仮定なしの

  `contHomCard_absGal_le : N_n(Γ_K) ≤ n · gcd(n, q−1)`

を得る。

## 帰着の骨格

`restrictInertia` の核は「惰性群 `I_K` の上で自明な連続指標」の全体である。
`π : Γ_K ↠ Gal(K^ur/K)`(`restrictNormalHom`)の核はちょうど `I_K`
(`IntermediateField.restrictNormalHom_ker`)なので、そのような `f` は
`f = g ∘ π` と**一意に**分解する(`MonoidHom.liftOfSurjective`)。あとは

* `g` が連続(核が開)であること、
* `f ↦ g` が単射であること(`π` が全射なので明らか)

を言えば `#ker(restrictInertia) ≤ #Hom_cont(Gal(K^ur/K), ℤ/n) = n` が出る。
★**単射で十分**で、全単射(=`Hom_cont` の同型)は要らない。上界だけが目的だから。

## `g` の連続性 —— 本ファイルの実質

`f` の核 `H` は開なので

1. `H` は閉(`Subgroup.isClosed_of_isOpen`)。よって無限次 Galois 対応
   (`InfiniteGalois.fixingSubgroup_fixedField`)で `H = (fixedField H).fixingSubgroup`。
2. `I_K ≤ H` から `fixedField H ≤ fixedField I_K = K^ur`
   (`InfiniteGalois.fixedField_fixingSubgroup`)。すなわち **`H` の固定体は不分岐**。
3. `fixedField H` は単項(在庫 `exists_adjoin_eq_fixedField`)なので、
   ある `y ∈ K^ur` について `H = K⟮y⟯.fixingSubgroup`。
4. `K^ur` の中の `K⟮y⟯` の固定部分群は `Gal(K^ur/K)` の**開**部分群であり
   (`IntermediateField.fixingSubgroup_isOpen`)、`π` で持ち上げると `H` に入るので
   `ker g` に含まれる。よって `ker g` は開(`Subgroup.isOpen_mono`)。

★2 が要である。「開部分群が惰性群を含む ⇔ その固定体が不分岐」という、
不分岐拡大の定義そのものを使わない言い換えになっている
——`IsUnramifiedAt` も `inertia` も経由しない。

## ★惰性群を `inertia` で書かない(`InertiaKummer.lean` の設計を引き継ぐ)

`Found/PGC/InertiaIdentification.lean` の `inertia` は `Interface` の
`SubgroupCorrespondence` / `ResidueCardinality` を経由し、経路 C の出口である
Corollary 1.3 と循環する。本ファイルも惰性群を常に
`(unramifiedClosure K).fixingSubgroup` と直接書き、`inertia` を参照しない。

## 逸脱の記録

* 惰性群を `inertia` ではなく `(unramifiedClosure K).fixingSubgroup` と書く(上記)。
* 「連続」を `IsOpen (ker f)` で表す(`Found/PGC/ContinuousHomCount.lean` の規約。
  `ZMod n` に `TopologicalSpace` インスタンスが無いための措置)。
* 原典(および Serre の局所類体論)は (C-q) を相互律から出すが、本経路は経由しない。
  これは `ResearchPaper/pgc-goal.md` に記録済みの逸脱である。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued

variable {p : ℕ} [Fact p.Prime]

/-! ## `Γ_K ↠ Gal(K^ur/K)` -/

/-- `K^ur` への制限射 `π : Γ_K → Gal(K^ur/K)`。

`Found/PGC/InertiaIdentification.lean` の `absGalQuotKerEquivUnramifiedGal` と同じ射だが、
そちらは商を作るために `inertia`(= `Interface` 経由)を通る。ここでは循環を避けるため
`restrictNormalHom` を直接名前で持つ。 -/
noncomputable abbrev restrictUnramified (K : PAdicLocalField p) :
    K.absGal →* (↥(unramifiedClosure K) ≃ₐ[K.carrier] ↥(unramifiedClosure K)) :=
  AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure) ↥(unramifiedClosure K)

/-- `π` は全射(`K^ur/K` は normal、`K̄/K` も normal)。 -/
theorem surjective_restrictUnramified (K : PAdicLocalField p) :
    Function.Surjective (restrictUnramified K) := by
  haveI := IsAlgClosure.normal K.carrier K.closure
  exact AlgEquiv.restrictNormalHom_surjective (F := K.carrier)
    (K₁ := ↥(unramifiedClosure K)) (E := K.closure)

/-- `π` の核はちょうど惰性群 `I_K = Gal(K̄/K^ur)`。 -/
theorem ker_restrictUnramified (K : PAdicLocalField p) :
    (restrictUnramified K).ker = (unramifiedClosure K).fixingSubgroup :=
  IntermediateField.restrictNormalHom_ker _

/-! ## 惰性群を含む開部分群は、`K^ur` の元 1 つの固定部分群 -/

/-- **★★★★★★★★惰性群を含む `Γ_K` の開部分群は `K^ur` の中で単項**。

`H` が開で `I_K ≤ H` なら、ある `y ∈ K^ur` について `H = K⟮y⟯.fixingSubgroup`。

証明は無限次 Galois 対応そのもの:`H` は開なので閉、よって
`H = (fixedField H).fixingSubgroup`(`InfiniteGalois.fixingSubgroup_fixedField`)。
`I_K ≤ H` は `fixedField H ≤ fixedField I_K = K^ur` を与える
(`InfiniteGalois.fixedField_fixingSubgroup`)。単項性は在庫
`exists_adjoin_eq_fixedField`(開部分群の固定体は `K` 上有限次で、原始元を持つ)。

退化の自己検査。

* `hH`(開)を落とすと**偽**——閉でない部分群の固定体は一般に有限次でなく、
  `fixingSubgroup (fixedField H) = H` も成り立たない。
* `hI`(惰性群を含む)を落とすと `y ∈ K^ur` が言えない。`H` が完全分岐拡大の
  固定部分群なら固定体は `K^ur` と `K` でしか交わらない。
* `H = ⊤` では `y ∈ K` が取れて両辺とも `⊤`(自明に真)。 -/
theorem exists_mem_unramifiedClosure_fixingSubgroup_eq (K : PAdicLocalField p)
    {H : Subgroup K.absGal} (hH : IsOpen (H : Set K.absGal))
    (hI : (unramifiedClosure K).fixingSubgroup ≤ H) :
    ∃ y : ↥(unramifiedClosure K),
      (IntermediateField.adjoin K.carrier ({(y : K.closure)} : Set K.closure)).fixingSubgroup
        = H := by
  haveI := isGalois_closure K
  obtain ⟨x, hx⟩ := exists_adjoin_eq_fixedField K H hH
  have hxfix : x ∈ IntermediateField.fixedField H := by
    rw [← hx]; exact IntermediateField.mem_adjoin_simple_self _ _
  have hxmem : x ∈ unramifiedClosure K := by
    have h2 : x ∈ IntermediateField.fixedField ((unramifiedClosure K).fixingSubgroup) := by
      rw [IntermediateField.mem_fixedField_iff] at hxfix ⊢
      exact fun σ hσ => hxfix σ (hI hσ)
    rwa [InfiniteGalois.fixedField_fixingSubgroup] at h2
  refine ⟨⟨x, hxmem⟩, ?_⟩
  show (IntermediateField.adjoin K.carrier ({x} : Set K.closure)).fixingSubgroup = H
  rw [hx]
  exact InfiniteGalois.fixingSubgroup_fixedField ⟨H, Subgroup.isClosed_of_isOpen H hH⟩

/-- **★★★★★★★★★★★★`π` を経由する連続指標の相棒も連続**。

`f : Γ_K →* ℤ/n` が連続(核が開)で惰性群の上で自明なら、`f = g ∘ π` をみたす
`g : Gal(K^ur/K) →* ℤ/n` も連続。

`ker f` を上の補題で `K⟮y⟯.fixingSubgroup`(`y ∈ K^ur`)と書き、`K^ur` の中の
`K⟮y⟯` の固定部分群を `ker g` の中に見つける。`τ` が `y` を固定すれば、`π` の
任意の持ち上げ `σ` も `y` を固定する(`restrictNormalHom_apply`)ので
`σ ∈ ker f`、すなわち `g τ = f σ = 1`。

退化の自己検査。

* `hI`(惰性群の上で自明)を落とすと `f` が `g ∘ π` の形にならないので、
  主張の前提 `hgf` 自体が置けない。
* `hf`(`f` が連続)を落とすと**偽**:`Gal(K^ur/K) ≅ Ẑ` の抽象的な指標は
  `n > 1` なら非可算個あり、そのどれも `π` と合成すれば `hgf` をみたす。 -/
theorem isOpen_ker_of_factors (K : PAdicLocalField p) {n : ℕ}
    {g : (↥(unramifiedClosure K) ≃ₐ[K.carrier] ↥(unramifiedClosure K)) →* Multiplicative (ZMod n)}
    {f : K.absGal →* Multiplicative (ZMod n)}
    (hf : IsOpen ((MonoidHom.ker f : Subgroup K.absGal) : Set K.absGal))
    (hI : (unramifiedClosure K).fixingSubgroup ≤ MonoidHom.ker f)
    (hgf : ∀ σ : K.absGal, g (restrictUnramified K σ) = f σ) :
    IsOpen ((MonoidHom.ker g : Subgroup (↥(unramifiedClosure K) ≃ₐ[K.carrier]
      ↥(unramifiedClosure K))) : Set (↥(unramifiedClosure K) ≃ₐ[K.carrier]
      ↥(unramifiedClosure K))) := by
  obtain ⟨y, hy⟩ := exists_mem_unramifiedClosure_fixingSubgroup_eq K hf hI
  haveI : FiniteDimensional K.carrier
      (IntermediateField.adjoin K.carrier ({y} : Set ↥(unramifiedClosure K))) :=
    IntermediateField.finiteDimensional_adjoin (fun z _ => Algebra.IsIntegral.isIntegral z)
  haveI : FiniteDimensional K.carrier
      (IntermediateField.adjoin K.carrier ({(y : K.closure)} : Set K.closure)) :=
    IntermediateField.finiteDimensional_adjoin (fun z _ => Algebra.IsIntegral.isIntegral z)
  have hle : (IntermediateField.adjoin K.carrier
      ({y} : Set ↥(unramifiedClosure K))).fixingSubgroup ≤ MonoidHom.ker g := by
    intro τ hτ
    obtain ⟨σ, rfl⟩ := surjective_restrictUnramified K τ
    rw [MonoidHom.mem_ker, hgf σ]
    have hτy : (restrictUnramified K σ) y = y :=
      (IntermediateField.mem_fixingSubgroup_iff _ _).mp hτ y
        (IntermediateField.mem_adjoin_simple_self _ _)
    have hσy : σ (y : K.closure) = (y : K.closure) := by
      rw [← AlgEquiv.restrictNormalHom_apply (unramifiedClosure K) σ y]
      exact congrArg Subtype.val hτy
    have hmem := mem_fixingSubgroup_adjoin_of_eq_self K (y : K.closure) σ hσy
    rw [hy] at hmem
    exact hmem
  exact Subgroup.isOpen_mono hle (IntermediateField.fixingSubgroup_isOpen _)

/-! ## `#ker(restrictInertia) ≤ n` -/

/-- `restrictInertia` の核に入る指標は、惰性群の上で値 `1`(核の定義の言い換え)。 -/
theorem eq_one_of_mem_ker_restrictInertia (K : PAdicLocalField p) {n : ℕ}
    {f : ↥(contHom K.absGal (Multiplicative (ZMod n)))} (hf : f ∈ (restrictInertia K n).ker)
    {σ : K.absGal} (hσ : σ ∈ (unramifiedClosure K).fixingSubgroup) :
    (f : K.absGal →* Multiplicative (ZMod n)) σ = 1 := by
  have h1 := MonoidHom.mem_ker.mp hf
  have h3 : ((restrictInertia K n) f :
      ↥((unramifiedClosure K).fixingSubgroup) →* Multiplicative (ZMod n)) ⟨σ, hσ⟩ = 1 := by
    rw [h1]; rfl
  exact h3

/-- **★★★★★★★★★★★★★★★★(F4)`#ker(restrictInertia) ≤ n`**。

`Found/PGC/InertiaKummer.lean::contHomCard_absGal_le_of_card_ker` が仮定に残していた
一点。惰性群の上で自明な連続指標 `f` は `f = g ∘ π` と分解し、`g` も連続
(`isOpen_ker_of_factors`)なので、`f ↦ g` は
`Hom_cont(Gal(K^ur/K), ℤ/n)`(F3 でちょうど `n` 個)への単射を与える。

退化の自己検査。

* `hn : n ≠ 0` を落とすと `ZMod 0 = ℤ` で右辺が `0` になり**偽**
  (左辺は自明指標のぶん少なくとも `1`)。
* 不等号を等号に強めることもできる(`π` が全射なので `f ↦ g` は全単射)が、
  上界しか要らないので単射で止めている。 -/
theorem card_ker_restrictInertia_le (K : PAdicLocalField p) {n : ℕ} (hn : n ≠ 0) :
    Nat.card ((restrictInertia K n).ker) ≤ n := by
  have hkerle : ∀ f : ↥((restrictInertia K n).ker),
      (restrictUnramified K).ker ≤ MonoidHom.ker (f.1 : K.absGal →* Multiplicative (ZMod n)) := by
    intro f σ hσ
    rw [ker_restrictUnramified] at hσ
    exact MonoidHom.mem_ker.mpr (eq_one_of_mem_ker_restrictInertia K f.2 hσ)
  have hgf : ∀ (f : ↥((restrictInertia K n).ker)) (σ : K.absGal),
      ((restrictUnramified K).liftOfSurjective (surjective_restrictUnramified K)
        ⟨(f.1 : K.absGal →* Multiplicative (ZMod n)), hkerle f⟩) (restrictUnramified K σ)
        = (f.1 : K.absGal →* Multiplicative (ZMod n)) σ :=
    fun f σ => MonoidHom.liftOfRightInverse_comp_apply _ _ _ _ σ
  have hopen : ∀ f : ↥((restrictInertia K n).ker),
      ((restrictUnramified K).liftOfSurjective (surjective_restrictUnramified K)
        ⟨(f.1 : K.absGal →* Multiplicative (ZMod n)), hkerle f⟩)
        ∈ contHom (↥(unramifiedClosure K) ≃ₐ[K.carrier] ↥(unramifiedClosure K))
          (Multiplicative (ZMod n)) := by
    intro f
    rw [mem_contHom]
    refine isOpen_ker_of_factors K f.1.2 (fun σ hσ => hkerle f ?_) (hgf f)
    rw [ker_restrictUnramified]
    exact hσ
  haveI : Finite ↥(contHom (↥(unramifiedClosure K) ≃ₐ[K.carrier] ↥(unramifiedClosure K))
      (Multiplicative (ZMod n))) :=
    Nat.finite_of_card_ne_zero (by
      rw [show Nat.card ↥(contHom (↥(unramifiedClosure K) ≃ₐ[K.carrier] ↥(unramifiedClosure K))
        (Multiplicative (ZMod n))) = contHomCard
        (↥(unramifiedClosure K) ≃ₐ[K.carrier] ↥(unramifiedClosure K)) n from rfl,
        contHomCard_unramifiedClosureGal K hn]
      exact hn)
  have hinj : Function.Injective
      (fun f : ↥((restrictInertia K n).ker) => (⟨_, hopen f⟩ :
        ↥(contHom (↥(unramifiedClosure K) ≃ₐ[K.carrier] ↥(unramifiedClosure K))
          (Multiplicative (ZMod n))))) := by
    intro a b hab
    have h1 : ((restrictUnramified K).liftOfSurjective (surjective_restrictUnramified K)
          ⟨(a.1 : K.absGal →* Multiplicative (ZMod n)), hkerle a⟩)
        = ((restrictUnramified K).liftOfSurjective (surjective_restrictUnramified K)
          ⟨(b.1 : K.absGal →* Multiplicative (ZMod n)), hkerle b⟩) :=
      congrArg Subtype.val hab
    refine Subtype.ext (Subtype.ext (MonoidHom.ext fun σ => ?_))
    rw [← hgf a σ, ← hgf b σ, h1]
  calc Nat.card ((restrictInertia K n).ker)
      ≤ Nat.card ↥(contHom (↥(unramifiedClosure K) ≃ₐ[K.carrier] ↥(unramifiedClosure K))
        (Multiplicative (ZMod n))) := Nat.card_le_card_of_injective _ hinj
    _ = n := contHomCard_unramifiedClosureGal K hn

/-! ## (C-q) の上界 -/

/-- **★★★★★★★★★★★★★★★★★★★★(G2)`N_n(Γ_K) ≤ n · gcd(n, q−1)`**(`p ∤ n`)。

経路 C(`ResearchPaper/pgc-goal.md`)の (C-q) の上界。`InertiaKummer.lean` の
`contHomCard_absGal_le_of_card_ker` から**仮定を外した形**である
(仮定つきの版は消していない)。

内訳は第一同型定理:

* 像の側 `≤ gcd(n, q−1)` —— `card_range_restrictInertia_le`(第 1033、F1+F2)。
* 核の側 `≤ n` —— `card_ker_restrictInertia_le`(本ファイル、F4)+
  `contHomCard_unramifiedClosureGal`(第 1035、F3)。

退化の自己検査。

* `hn : ¬ p ∣ n` を落とすと**偽**。`n = p` では惰性群からの連続指標が
  野性的惰性群のぶんだけ増え、(F1) の「位数 `n` の巡回群」が崩れる。
* `n ∣ q − 1` のとき右辺は `n²` で、下界(在庫 `contHomCard_absGal_of_dvd`)と
  一致する。すなわちこの上界はその範囲で**最良**である。
* `n = 1` では両辺 `1`(自明に真)。 -/
theorem contHomCard_absGal_le (K : PAdicLocalField p) {n : ℕ} (hn : ¬ p ∣ n) :
    contHomCard K.absGal n ≤ n * Nat.gcd n (Nat.card 𝓀[K.carrier] - 1) :=
  contHomCard_absGal_le_of_card_ker K hn
    (card_ker_restrictInertia_le K (by rintro rfl; exact hn (dvd_zero p)))

end ABC3.Found.PGC
