import ABC3.Found.PGC.ArtinEquivarianceProof

/-!
# `σ_g` が Lubin-Tate 塔を運ぶ段 —— (i) の完成と (ii) の材料

直前の節点(`ArtinEquivarianceProof.lean`)は、[pGC] Proposition 1.1 に残った壁を
`ArtinUnitEquivariance`(Artin 写像の**単数部分**の `σ`-同変性)まで落とし、
その申し送りとして 3 段を名指しした:

> **(i)** §13(捩れ点の輸送)を `IntermediateField.adjoin` に持ち上げて
>   `σ_g (K_{f,n}) = K_{f^{σ_g},n}`
> **(ii)** `ρ_{f^{σ_g}}(g τ g^{-1}) = σ_g(ρ_f(τ))`
> **(iii)** 素元非依存性で `π` に戻す

本ファイルは **(i) を完成させ**、**(ii) の材料を 3 本用意した**。
★★**`ArtinUnitEquivariance` は出ていない。** 出ていないものを出たと書かない。

## ★★何がこのファイルで出たか(★測定結果)

| 段 | 宣言 | 状態 |
|---|---|---|
| (i-a) | `image_iteratedLubinTateTorsionPoints` | ★**証明した**(`Φ '' Λ_{f,n} = Λ_{f^φ,n}`、包含ではなく**等号**) |
| (i-b) | `image_lubinTateTorsionSet` | ★**証明した**(`Φ '' Λ_∞ = Λ'_∞`) |
| (i-c) | `image_lubinTateLevelField` | ★★**証明した** —— これが (i) 本体 |
| (i-d) | `image_lubinTateClosure` | ★**証明した**(塔全体の版) |
| (i-e) | `image_lubinTateLevelField_fixedFieldClosureAut` | ★**証明した**(`Γ_F` の `g` に代入した形) |
| (ii-a) | `map_LubinTateAction` | ★★**証明した** —— `([a]_f)^φ = [φ a]_{f^φ}`(冪級数の層) |
| (ii-b) | `aeval_powerSeries_comm_twist` | ★★**証明した** —— 半線型な連続環準同型は `PowerSeries.aeval` を運ぶ |
| (ii-c) | `eq_of_conj_spec_of_existsUnique` | ★**証明した** —— 一意性から `ρ'(Φ τ Φ^{-1}) = φ(ρ τ)` を出す抽象核 |
| (橋) | `absGalConjCME_apply` | ★★**証明した(`rfl`)** —— `Ψ_g` は `Φ_g` による共役**そのもの** |
| (ii) | `ρ_{f^σ}(Φ τ Φ^{-1}) = σ(ρ_f τ)` | ★★**出ていない**(下記 3 点が足りない) |
| (iii) | 素元非依存(写像の等式) | ★★**手を付けていない** |

## ★★★残っているものを名指しで(★次の節点はここから持てる)

(ii) を閉じるのに足りないのは **3 つだけ**で、どれも本ファイルの
`aeval_powerSeries_comm_twist` に**代入する**段である:

1. ★**`Φ` を `adjoinIntegers K x → adjoinIntegers K (Φ x)` に制限する段**。
   `Φ '' K⟮x⟯ = K⟮Φ x⟯` は本ファイルの `image_adjoin_of_semilinear_equiv` に
   `S := {x}` を入れれば出る。残るのは
   (a) `Φ` がノルムを保つこと(`Φ` は `ℚ_p`-代数同型なのでスペクトルノルムを保つ。
       木の `norm_algEquiv_eq` は `≃ₐ[K.carrier]` 専用で、★**半線型版は無い**)、
   (b) 制限が連続であること(`LubinTateActionEquivariance.lean` の
       `continuous_algEquivRestrictSelf` と同じ形。★ただしあちらは
       `σ(x) ∈ K⟮x⟯` に留まる場合専用で、半線型では**留まらない**
       ——`Φ x` は `f^φ` の捩れ点であって `f` の捩れ点ではない)。
   ★★**これが木が「cross-point instance bridging」と呼んで避けてきた段である。**
   半線型では避けられない(避けると `φ = id` になって主張が空になる)。
2. ★**`reciprocityUnits` の「捩れ点の上での記述」**。木にあるのは
   `mem_principalUnits_reciprocityUnits_iff`(核の記述)だけで、
   ★**`σ (x_m) = [reciprocityUnits σ]_f (x_m)` という形の spec は無い**
   (`reciprocityMap_spec` を `unitsEquivCompatibleUnits` /
   `principalUnitsQuotientEquiv` を通して降ろす段が要る)。
3. ★**生成元列の付け替え**。`reciprocityMapLimit` は `psiGenSeq`(選択された
   compatible 生成元列)で定義されているが、`f^φ` 側の `psiGenSeq` は
   `Φ '' (f 側の psiGenSeq)` とは限らない。★**一意性
   (`existsUnique_unitActionQuotient_eq_algEquiv`)で任意の生成元に移す段**が要る。

(iii)(素元非依存)は本ファイルでは**手を付けていない**。★**必要である**ことは
測って確かめた: `σ_g(π)` は `π` と同伴な**別の**素元であり、`ρ_f` と `ρ_{f^{σ_g}}`
は `Γ_K` の上では一致しない(`Art(π)` の像が違う)。一致するのは**慣性部分**
(`Ẑ` 成分が `1` の部分、すなわち `p^n` 捩れが乗る部分)だけで、それが
原典 Corollary 4.9 の `j = 0` の場合である。★**`ArtinUnitEquivariance` が
`p^n` 捩れの上でしか要求していないことが、ここで効く。**

## ★★測ってひっくり返ったこと(★段取りとの差分)

* ★★**`Ψ_g`(`absGalConjCME`)は `Φ_g := Θ^{-1} ∘ g ∘ Θ` による共役そのもので、
  等式は `rfl` で通った**(`absGalConjCME_apply`)。
  申し送りは「#59 の (a) か (e) になる見込み」と読んでいたが、
  ★**#59 には 1 度も当たらなかった** —— `closureEquivFixedField` を
  `RingEquiv` として 1 本に潰した時点で中間体の層が消える(#59 の定型 (d))。
* ★★**(i) は「濃度で押す」と 1 行で閉じた。** 直前の波の
  `map_mem_iteratedLubinTateTorsionPoints` は**包含**しか与えないが、
  `Φ` は単射で `|Λ_{f,n}| = |Λ_{f^φ,n}| = q^n`(`card_iteratedLubinTateTorsionPoints`)
  なので `Finset.eq_of_subset_of_card_le` で**等号**になる。
  ★逆向きの写像(`φ^{-1}` で戻す)を作る必要が無かった。
* ★★**`adjoin` への持ち上げに `IntermediateField.adjoin_map` は使えない**
  (あれは基礎体を固定する `AlgHom` 専用で、`Φ` は `K.carrier` 上線型でない)。
  ★代わりに `Subfield.closure` に落とすと 3 行で出る
  (`subfield_map_closure` + `IntermediateField.adjoin_toSubfield`)。
  ★**`Subfield.map_closure` は mathlib に無い**
  (★どう測ったか: `grep -n "Subfield.map" .cache/mathlib-index.txt | head -20`。
  `map_bot` / `map_comap_eq` / `map_comap_eq_self` / `map_comap_eq_self_of_surjective` /
  `map_iInf` / `map_iSup` / `map_inf` / `map_le_iff_le_comap` / `map_map` /
  `map_mem_map` / `map_sup` が出るが `map_closure` は出ない。
  ★`Subfield.closure_le` 1 本で両向きが出るので自分で書いた)。
* ★**`LubinTateEndoTwisted` は今回も消費していない。** 本ファイルの
  `map_LubinTateAction`(`[a]_f` の係数ひねり)は在庫
  `powerSeries_uniqueness`(`LubinTateUniqueness.lean`、★**untwisted 版**)と
  `map_subst_powerSeries`(`DworkThetaStep2.lean`)の 2 本で出た。
  ★★**ねじれ版の一意性(`powerSeries_uniqueness_twisted`)は要らなかった** ——
  `f` と `f^φ` は**同じ環 `𝒪_K` の上**にあり、`ϕ = id` で足りるからである。
  ★これは 3 波連続で「`LubinTateEndoTwisted` を消費するはず」という見立てが
  外れたことになる。次の節点は**見立てを立てずに型で引く**方がよい。
* ★**素元非依存性も今回は消費していない**(上記 (iii))。

## ★退化の自己検査

1. ★**`Ẑ` 成分は本ファイルにも現れない。** 扱っているのは `𝒪_K^×` 側
   (= 捩れ点・Lubin-Tate 塔)だけである。直前の波の測定を壊していない。
2. ★**`ArtinUnitEquivariance` の主張も `artinUnitEquivariance_iff_cyclotomeConj`
   も触っていない**(本ファイルは `ArtinEquivarianceProof` を import するだけで、
   その中の宣言を 1 つも書き換えていない)。同値性は壊れていない。
3. ★**易しい半分(`g ∈ S`)は独立に成り立ったままである。** 本ファイルの
   `fixedFieldClosureAut` は `g ∈ S` のとき恒等とは限らない(`Θ^{-1} g Θ` は
   `K.closure` の上では自明でない)が、`σ_g` の側は `fixedFieldAut_of_mem` で
   恒等になる。★したがって `image_lubinTateLevelField` は `g ∈ S` のとき
   「`Φ_g` は `K(Λ_{f,n})` を**自分自身**に写す」を言う(`f^{σ_g} = f`)。
   ★これは既知の事実(`K(Λ_{f,n})/K` は正規)と整合する。
4. ★**`Φ` に `K`-線型性を仮定していない。** 仮定すると `φ = id` になり
   (`σ` が `K.carrier` の恒等になり)主張が空になる。★本ファイルの
   すべての抽象核でこの点を明示した。
5. ★**`n = 0` でも成り立つ**(`Λ_0 = {0}`、両辺とも `K`)。
6. ★**`ℕ∞` の切り詰め引き算・除算を 1 つも書いていない**(#102)。

## 逸脱の記録

* 新しい仮定を 1 つも置いていない。本ファイルの定理はすべて
  `sorry` 無し・仮説付きの**無条件**な主張である。
* 原典 Lemma 4.6 は `θ ∈ Θ^{L,×}_{π,π′}` による `[θ]` の輸送を扱うが、
  本ファイルが扱うのは **`σ`(体の自己同型)による輸送**である。
  ★どちらも「`f` を取り替えると塔が動く」という同じ現象の 2 つの面だが、
  ★**同じ主張ではない**(`[θ]` は `𝒪`-線型、`σ` は半線型)。
  `.src` を Lemma 4.6 に付けたのは「体の等式 `L^m_f = L^m_{f′}` を出す段」
  という**役割**が同じだからで、証明は原典と別である(原典は `[θ]` と
  `[θ^{-1}]` の両方が `𝒪_L[[X]]` にあることを使うが、本ファイルは
  捩れ点の**濃度**を使う)。★この逸脱をここに記録する。
* `fixedFieldClosureAut` は `closureEquivFixedField`(`IsAlgClosure.equiv`、
  選択公理)を経由する。★同型の**取り方**には依存するが、
  `absGalConjCME_apply` が示すとおり、`Ψ_g` との一致は取り方に依らず `rfl` で出る。

## ★測定の記録(★次の agent が再実行できる形で)

* ★**mathlib を先に引いた**:
  `grep -n "Subfield.map" .cache/mathlib-index.txt | head -20`
  (`map_closure` は出ない ——上の節に出た一覧が全部である)、
  `grep -nE "	IntermediateField\.(adjoin_toSubfield|adjoin_map|restrictScalars_adjoin)	"`
  → `IntermediateField.adjoin_toSubfield` を引き当て(これが (i) の要)、
  `grep -nE "	PowerSeries\.(trunc_map)"` → `PowerSeries.trunc_map` を引き当て
  (`aeval_powerSeries_comm_twist` の truncation 段)。
  ★**「mathlib に無い」と書いたのは `Subfield.map_closure` の 1 件だけ**で、
  それは上のコマンドで**測った**結果である。
* ★**#158(同名で書いて `already been declared`)**: `subfield_map_closure` /
  `map_LubinTateAction` / `aeval_powerSeries_comm_twist` /
  `fixedFieldClosureAut` の 4 つを試して **0 件**。
  一方 `grep -n "powerSeries_uniqueness" .cache/decl-index.txt` は
  `powerSeries_uniqueness`(`LubinTateUniqueness.lean:174`)を**引き当てた**
  ——★`map_LubinTateAction` の一意性はこれ 1 本で済んだ(自分で書いていない)。
  同じく `grep -n "LubinTateAction" .cache/decl-index.txt | grep -i "map"` は
  12 件返すが**どれも「係数写像との可換性」ではない**(`reciprocityMap` 系と
  `_map_residue` 系)。★「`[a]_f` の係数ひねりの自然性」が木に無いことを
  確かめてから書いた。
* ★**#59 には当たっていない。** 中間体は `K.closure` の中の `IntermediateField`
  1 層だけで、`restrictNormalHom` も `restrictScalars`(体の)も 1 度も書いていない。
  `closureEquivFixedField` は `RingEquiv` に潰してある(定型 (d))。
* ★**往復回数**: `leanfile.mjs` で 17 往復(import だけの基準 9.6 秒、
  各回 9.6〜11.5 秒 ⇒ 周辺は 0.0〜1.9 秒)。
  ★**12 個の断片のうち 9 個が一発で通った**(直らなかった 3 個はいずれも
  coercion か暗黙引数の問題で、数学ではない ——#205 / #206 に記録した)。
  ★★**共有 REPL(`lean_start`)は別の import で走っていた(`lean_status` で確認)ので
  1 度も触らず、最初から `leanfile.mjs` に寄せた。**
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC

open scoped NormedField Valued Classical IsMulCommutative

open scoped PowerSeries.WithPiTopology

/-! ## 1. 抽象核 A —— 部分体の生成・添加体と半線型写像

★この節に体論以外の語彙(付値・分岐・Galois・Lubin-Tate)は 1 つも出てこない。 -/

section AbstractCoreA

/-- ★**抽象核 A1** —— 部分体の生成は環準同型の像と可換。

★mathlib に `Subfield.map_closure` は無い(`Subfield.map_sup` / `map_iSup` /
`map_inf` / `map_iInf` / `map_bot` / `map_comap_eq` / `map_le_iff_le_comap` /
`map_map` / `map_mem_map` はある)。両向きとも `Subfield.closure_le` 1 本で出る。 -/
theorem subfield_map_closure {K L : Type*} [Field K] [Field L] (f : K →+* L) (T : Set K) :
    (Subfield.closure T).map f = Subfield.closure (f '' T) := by
  refine le_antisymm ?_ (Subfield.closure_le.mpr ?_)
  · rw [Subfield.map_le_iff_le_comap]
    exact Subfield.closure_le.mpr fun t ht => Subfield.subset_closure ⟨t, ht, rfl⟩
  · rintro _ ⟨t, ht, rfl⟩
    exact ⟨t, Subfield.subset_closure ht, rfl⟩

/-- ★★**抽象核 A2** —— 基礎体を(集合として)保つ半線型な `Φ` は
`k(S)` を `k(Φ '' S)` に写す。

★★`Φ` に `k`-線型性を仮定していないことが要点である。仮定すると
`IntermediateField.adjoin_map` がそのまま使えるが、そのときは
本ファイルの主張が空になる(`σ = id` の場合しか扱えない)。
★代わりに `IntermediateField.adjoin_toSubfield` で `Subfield.closure` に落とす。 -/
theorem image_adjoin_of_semilinear {k Ω : Type*} [Field k] [Field Ω] [Algebra k Ω]
    (Φ : Ω →+* Ω) (hk : Φ '' Set.range (algebraMap k Ω) = Set.range (algebraMap k Ω))
    (S : Set Ω) :
    Φ '' (IntermediateField.adjoin k S : Set Ω)
      = (IntermediateField.adjoin k (Φ '' S) : Set Ω) := by
  have h1 : ((IntermediateField.adjoin k S).toSubfield.map Φ : Set Ω)
      = (Subfield.closure (Φ '' (Set.range (algebraMap k Ω) ∪ S)) : Set Ω) := by
    rw [IntermediateField.adjoin_toSubfield, subfield_map_closure]
  rw [Set.image_union, hk] at h1
  rw [show ((IntermediateField.adjoin k (Φ '' S) : IntermediateField k Ω) : Set Ω)
      = ((IntermediateField.adjoin k (Φ '' S)).toSubfield : Set Ω) from rfl,
    IntermediateField.adjoin_toSubfield, ← h1]
  rfl

/-- ★**抽象核 A3** —— `Φ` が `σ` について半線型なら基礎体の像は基礎体。 -/
theorem image_range_algebraMap_of_semilinear {k Ω : Type*} [Field k] [Field Ω] [Algebra k Ω]
    (Φ : Ω →+* Ω) (σ : k ≃+* k)
    (hΦ : ∀ a : k, Φ (algebraMap k Ω a) = algebraMap k Ω (σ a)) :
    Φ '' Set.range (algebraMap k Ω) = Set.range (algebraMap k Ω) := by
  refine Set.Subset.antisymm ?_ ?_
  · rintro _ ⟨_, ⟨a, rfl⟩, rfl⟩
    exact ⟨σ a, (hΦ a).symm⟩
  · rintro _ ⟨b, rfl⟩
    exact ⟨algebraMap k Ω (σ.symm b), ⟨σ.symm b, rfl⟩, by rw [hΦ, σ.apply_symm_apply]⟩

/-- ★**抽象核 A2 の `RingEquiv` 版**(coercion を `⇑Φ` の形に揃える)。

★`⇑↑Φ`(`RingEquiv → RingHom → 関数`)と `⇑Φ` は defeq だが `rw` は
**構文的一致**しか見ないので、消費側の形に合わせた版を用意する(#200)。 -/
theorem image_adjoin_of_semilinear_equiv {k Ω : Type*} [Field k] [Field Ω] [Algebra k Ω]
    (Φ : Ω ≃+* Ω) (σ : k ≃+* k)
    (hΦ : ∀ a : k, Φ (algebraMap k Ω a) = algebraMap k Ω (σ a)) (S : Set Ω) :
    Φ '' (IntermediateField.adjoin k S : Set Ω)
      = (IntermediateField.adjoin k (Φ '' S) : Set Ω) :=
  image_adjoin_of_semilinear (Φ : Ω →+* Ω)
    (image_range_algebraMap_of_semilinear (Φ : Ω →+* Ω) σ hΦ) S

end AbstractCoreA

/-! ## 2. 抽象核 B —— 半線型輸送による「作用の記述」の運搬

★★原典 Proposition 4.7(ii) の

> If α ∈ µ×f,m, then σ(α) ∈ µ(j),×f,m , hence σ(α) = [xπj](α) for a unique x mod 1+pm

の、`σ` を**体の自己同型で共役する**版である。
★この節にも体・付値・分岐・Lubin-Tate の語彙は 1 つも出てこない。 -/

section AbstractCoreB

/-- ★★**抽象核 B1** —— 半線型な輸送 `Φ` は
「`τ` が `M` の上で `br a`」を「`Φ τ Φ^{-1}` が `Φ '' M` の上で `br' (φ a)`」に運ぶ。

★`Φ` が全単射である必要は無く、右逆 `Ψ` があればよい。
★`M` が部分加群である必要も、`br` が合成則を満たす必要も無い。 -/
theorem conj_spec_transfer {S S' N N' : Type*}
    (br : S → N → N) (br' : S' → N' → N') (φ : S → S')
    (Φ : N → N') (Ψ : N' → N) (hΨ : ∀ β : N', Φ (Ψ β) = β)
    (hsemi : ∀ (a : S) (α : N), Φ (br a α) = br' (φ a) (Φ α))
    (τ : N → N) (a : S) (M : Set N) (hτ : ∀ α ∈ M, τ α = br a α)
    (hMΨ : ∀ β ∈ Φ '' M, Ψ β ∈ M) :
    ∀ β ∈ Φ '' M, Φ (τ (Ψ β)) = br' (φ a) β := by
  intro β hβ
  rw [hτ _ (hMΨ β hβ), hsemi, hΨ]

/-- ★★★**抽象核 B2 —— (ii) の論理そのもの**:
`ρ` の記述の一意性から `ρ'(Φ τ Φ^{-1}) = φ(ρ τ)` が出る。

記号の対応(原典 Corollary 4.9 の証明の 1 行):

| 原典 | ここ |
|---|---|
| `µ_{f,m}` | `M` |
| `µ_{f^σ,m}` | `Φ '' M` |
| `[a]_f` / `[a]_{f^σ}` | `br` / `br'` |
| `σ`(係数への作用) | `φ` |
| `ρ_f(τ)` / `ρ_{f^σ}(Φ τ Φ^{-1})` | `rep q` / `rep' q'` |

★`huniq` は「`Φ '' M` の上で同じ作用をする係数は等しい」——
原典の "for a unique x mod 1+p^m" に当たる。 -/
theorem eq_of_conj_spec_of_existsUnique {S S' N N' Q Q' : Type*}
    (br : S → N → N) (br' : S' → N' → N') (φ : S → S')
    (Φ : N → N') (Ψ : N' → N) (hΨ : ∀ β : N', Φ (Ψ β) = β)
    (hsemi : ∀ (a : S) (α : N), Φ (br a α) = br' (φ a) (Φ α))
    (M : Set N) (hMΨ : ∀ β ∈ Φ '' M, Ψ β ∈ M)
    (rep : Q → S) (rep' : Q' → S') (τ : N → N)
    (q : Q) (hq : ∀ α ∈ M, τ α = br (rep q) α)
    (q' : Q') (hq' : ∀ β ∈ Φ '' M, Φ (τ (Ψ β)) = br' (rep' q') β)
    (huniq : ∀ s s' : S', (∀ β ∈ Φ '' M, br' s β = br' s' β) → s = s') :
    rep' q' = φ (rep q) :=
  huniq _ _ fun β hβ => by
    rw [← hq' β hβ, conj_spec_transfer br br' φ Φ Ψ hΨ hsemi τ (rep q) M hq hMΨ β hβ]

end AbstractCoreB

/-! ## 3. 抽象核 C —— 半線型な連続環準同型と `PowerSeries.aeval`

★★木の `PowerSeriesAevalComm.lean` は `σ : S₁ →ₐ[A] S₂`(**係数を固定**する
代数準同型)の場合しか持っていない。半線型な場合が要る
——これが「点で評価する層」の**作用側**である。 -/

section AbstractCoreC

/-- ★**抽象核 C0** —— `ArtinEquivarianceProof.lean` の `hom_eval₂_twist` の
「始域と終域が違う」版。 -/
theorem hom_eval₂_twist' {A B₁ B₂ : Type*} [CommRing A] [CommRing B₁] [CommRing B₂]
    (φ : A →+* A) (Φ : B₁ →+* B₂) (ι₁ : A →+* B₁) (ι₂ : A →+* B₂)
    (h : ∀ a, Φ (ι₁ a) = ι₂ (φ a)) (P : Polynomial A) (x : B₁) :
    Φ (Polynomial.eval₂ ι₁ x P) = Polynomial.eval₂ ι₂ (Φ x) (Polynomial.map φ P) := by
  rw [Polynomial.hom_eval₂, Polynomial.eval₂_map]
  congr 1
  exact RingHom.ext h

/-- ★★★★**抽象核 C —— 半線型な連続環準同型は `PowerSeries.aeval` を
「係数をひねった冪級数の評価」に運ぶ**:

  `Φ (g(z)) = (g^φ)(Φ z)`。

★分岐・付値・Galois・Lubin-Tate の語彙は 1 つも出てこない。純粋な位相環論である。
★★`Φ` に `A`-線型性を仮定していない(仮定すると `φ = id` になる)。
★木の `algHom_aeval_powerSeries_comm'`(`PowerSeriesAevalComm.lean`)は
`σ : S₁ →ₐ[A] S₂` 専用で、この形は無かった。
証明は truncation の極限で、多項式の段が抽象核 C0 になるだけである。 -/
theorem aeval_powerSeries_comm_twist {A S₁ S₂ : Type*} [CommRing A] [CommRing S₁] [CommRing S₂]
    [Algebra A S₁] [Algebra A S₂]
    [UniformSpace A] [IsUniformAddGroup A] [IsTopologicalSemiring A]
    [UniformSpace S₁] [IsUniformAddGroup S₁] [T2Space S₁] [CompleteSpace S₁] [IsTopologicalRing S₁]
    [IsLinearTopology S₁ S₁] [ContinuousSMul A S₁]
    [UniformSpace S₂] [IsUniformAddGroup S₂] [T2Space S₂] [CompleteSpace S₂] [IsTopologicalRing S₂]
    [IsLinearTopology S₂ S₂] [ContinuousSMul A S₂]
    (φ : A →+* A) (Φ : S₁ →+* S₂) (hΦcont : Continuous Φ)
    (hsemi : ∀ a : A, Φ (algebraMap A S₁ a) = algebraMap A S₂ (φ a))
    (g : PowerSeries A) {z : S₁} (hz : PowerSeries.HasEval z)
    (hz' : PowerSeries.HasEval (Φ z)) :
    Φ (PowerSeries.aeval hz g) = PowerSeries.aeval hz' (PowerSeries.map φ g) := by
  have htrunc : ∀ N : ℕ,
      Φ (PowerSeries.aeval hz ((PowerSeries.trunc N g : Polynomial A) : PowerSeries A))
        = PowerSeries.aeval hz'
            ((PowerSeries.trunc N (PowerSeries.map φ g) : Polynomial A) : PowerSeries A) := by
    intro N
    rw [PowerSeries.aeval_coe, PowerSeries.aeval_coe, PowerSeries.trunc_map,
      Polynomial.aeval_def, Polynomial.aeval_def]
    exact hom_eval₂_twist' φ Φ (algebraMap A S₁) (algebraMap A S₂) hsemi _ z
  have h1 : Filter.Tendsto
      (fun N => Φ (PowerSeries.aeval hz ((PowerSeries.trunc N g : Polynomial A) : PowerSeries A)))
      Filter.atTop (nhds (Φ (PowerSeries.aeval hz g))) :=
    hΦcont.continuousAt.tendsto.comp
      ((PowerSeries.continuous_aeval hz).continuousAt.tendsto.comp
        (PowerSeries.WithPiTopology.tendsto_trunc_atTop A g))
  have h2 : Filter.Tendsto
      (fun N => PowerSeries.aeval hz'
        ((PowerSeries.trunc N (PowerSeries.map φ g) : Polynomial A) : PowerSeries A))
      Filter.atTop (nhds (PowerSeries.aeval hz' (PowerSeries.map φ g))) :=
    (PowerSeries.continuous_aeval hz').continuousAt.tendsto.comp
      (PowerSeries.WithPiTopology.tendsto_trunc_atTop A (PowerSeries.map φ g))
  rw [funext htrunc] at h1
  exact tendsto_nhds_unique h1 h2

end AbstractCoreC

/-! ## 4. 冪級数層 —— `[a]_f` の係数ひねりに関する自然性 -/

def map_LubinTateAction.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 6, item := "Proposition 3.5", sectionId := "prop-3-5" }

/-- ★★★★**`([a]_f)^φ = [φ a]_{f^φ}`** —— Lubin-Tate の `𝒪`-作用は
係数のひねりと可換。

★★`f` と `f^φ` は**同じ環 `A` の上**にあるので、ねじれ版の一意性
(`LubinTateEndoTwisted.lean` の `powerSeries_uniqueness_twisted`)は要らない
——在庫 `powerSeries_uniqueness`(untwisted、`LubinTateUniqueness.lean`)に
`h := f^φ` を入れるだけである。関数等式は `map_subst_powerSeries`
(`DworkThetaStep2.lean`、係数写像は代入と可換)を 2 回使って運ぶ。 -/
theorem map_LubinTateAction {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    (φ : A ≃+* A)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (a : A) :
    PowerSeries.map (φ : A →+* A) (LubinTateAction hq hπmax f hf0 hf1 hf a)
      = LubinTateAction hq (maximalIdeal_eq_span_map φ hπmax)
          (PowerSeries.map (φ : A →+* A) f) (coeff_zero_map_twist _ hf0)
          (coeff_one_map_twist _ hf1) (map_residue_map_twist φ hf) (φ a) := by
  have hHSf : PowerSeries.HasSubst f := by
    show IsNilpotent (PowerSeries.constantCoeff f)
    rw [show PowerSeries.constantCoeff f = 0 by simpa using hf0]; exact IsNilpotent.zero
  have hHSa : PowerSeries.HasSubst (LubinTateAction hq hπmax f hf0 hf1 hf a) := by
    show IsNilpotent (PowerSeries.constantCoeff _)
    rw [constantCoeff_LubinTateAction]; exact IsNilpotent.zero
  refine powerSeries_uniqueness (maximalIdeal_eq_span_map φ hπmax)
    (fun h => hπne0 (φ.injective (by rw [h, map_zero]))) (h := PowerSeries.map (φ : A →+* A) f)
    (by rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]
        exact coeff_zero_map_twist (φ : A →+* A) hf0)
    (coeff_one_map_twist (φ : A →+* A) hf1) ?_ ?_ ?_ ?_ ?_
  · rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply, PowerSeries.coeff_map,
      PowerSeries.coeff_zero_eq_constantCoeff_apply, constantCoeff_LubinTateAction, map_zero]
  · exact constantCoeff_LubinTateAction _ _ _ _ _ _ _
  · rw [PowerSeries.coeff_map, coeff_one_LubinTateAction, coeff_one_LubinTateAction]
    rfl
  · rw [← map_subst_powerSeries (φ : A →+* A) hHSa, ← map_subst_powerSeries (φ : A →+* A) hHSf,
      LubinTateAction_functional_equation]
  · exact (LubinTateAction_functional_equation hq (maximalIdeal_eq_span_map φ hπmax)
      (PowerSeries.map (φ : A →+* A) f) (coeff_zero_map_twist _ hf0)
      (coeff_one_map_twist _ hf1) (map_residue_map_twist φ hf) (φ a)).symm

variable {p : ℕ} [Fact p.Prime]

/-! ## 5. 具体層 1 —— 半線型性を整数環へ落とす -/

/-- ★`Φ` が `σ`(体の同型)について半線型なら、整数環の上でも
`integerRingEquiv σ` について半線型。★在庫 `integerRingEquiv`
(`ResidueCardinality.lean`)に `IsScalarTower` を噛ませるだけ。 -/
theorem semilinear_integerRingEquiv (K : PAdicLocalField p) (σ : K.carrier ≃ₐ[ℚ_[p]] K.carrier)
    (Φ : K.closure →+* K.closure)
    (hΦ : ∀ a : K.carrier, Φ (algebraMap K.carrier K.closure a)
      = algebraMap K.carrier K.closure (σ a))
    (a : 𝒪[K.carrier]) :
    Φ (algebraMap 𝒪[K.carrier] K.closure a)
      = algebraMap 𝒪[K.carrier] K.closure (integerRingEquiv σ a) := by
  rw [IsScalarTower.algebraMap_apply 𝒪[K.carrier] K.carrier K.closure,
    IsScalarTower.algebraMap_apply 𝒪[K.carrier] K.carrier K.closure, hΦ]
  rfl

/-! ## 6. ★★★★具体層 2 —— (i) 捩れ点と Lubin-Tate 塔の輸送 -/

open scoped Classical in
/-- ★★**`Φ` は `Λ_{f,n}` を `Λ_{f^φ,n}` の上へ写す** —— 直前の波の
`map_mem_iteratedLubinTateTorsionPoints`(包含)を**等号**にした形。

★★等号は「単射 + 濃度」で出る:`|Λ_{f,n}| = |Λ_{f^φ,n}| = q^n`
(`card_iteratedLubinTateTorsionPoints`)。★逆向きの写像を作る必要は無い。 -/
theorem image_iteratedLubinTateTorsionPoints (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal 𝒪[K.carrier]) 𝒪[K.carrier]]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField 𝒪[K.carrier]) pp]
    [Fintype (IsLocalRing.ResidueField 𝒪[K.carrier])]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField 𝒪[K.carrier]) = pp ^ ff)
    (φ : 𝒪[K.carrier] ≃+* 𝒪[K.carrier])
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries 𝒪[K.carrier]) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (Φ : K.closure ≃+* K.closure)
    (hcompat : ∀ a : 𝒪[K.carrier], Φ (algebraMap 𝒪[K.carrier] K.closure a)
      = algebraMap 𝒪[K.carrier] K.closure (φ a)) :
    Finset.image Φ (iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
      = iteratedLubinTateTorsionPoints K hq (maximalIdeal_eq_span_map φ hπmax)
          (fun h => hπne0 (φ.injective (by rw [h, map_zero])))
          (PowerSeries.map (φ : 𝒪[K.carrier] →+* 𝒪[K.carrier]) f) (coeff_zero_map_twist _ hf0)
          (coeff_one_map_twist _ hf1) (map_residue_map_twist φ hf) n := by
  refine Finset.eq_of_subset_of_card_le ?_ (le_of_eq ?_)
  · intro y hy
    obtain ⟨x, hx, rfl⟩ := Finset.mem_image.mp hy
    exact map_mem_iteratedLubinTateTorsionPoints K hq φ hπmax hπne0 f hf0 hf1 hf n
      (Φ : K.closure →+* K.closure) hcompat hx
  · rw [Finset.card_image_of_injective _ Φ.injective,
      card_iteratedLubinTateTorsionPoints, card_iteratedLubinTateTorsionPoints]

/-- ★★`Φ '' Λ_∞ = Λ'_∞`(段ごとの像を合併するだけ)。 -/
theorem image_lubinTateTorsionSet (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal 𝒪[K.carrier]) 𝒪[K.carrier]]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField 𝒪[K.carrier]) pp]
    [Fintype (IsLocalRing.ResidueField 𝒪[K.carrier])]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField 𝒪[K.carrier]) = pp ^ ff)
    (φ : 𝒪[K.carrier] ≃+* 𝒪[K.carrier])
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries 𝒪[K.carrier]) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f = PowerSeries.X ^ (pp ^ ff))
    (Φ : K.closure ≃+* K.closure)
    (hcompat : ∀ a : 𝒪[K.carrier], Φ (algebraMap 𝒪[K.carrier] K.closure a)
      = algebraMap 𝒪[K.carrier] K.closure (φ a)) :
    Φ '' lubinTateTorsionSet K hq hπmax hπne0 f hf0 hf1 hf
      = lubinTateTorsionSet K hq (maximalIdeal_eq_span_map φ hπmax)
          (fun h => hπne0 (φ.injective (by rw [h, map_zero])))
          (PowerSeries.map (φ : 𝒪[K.carrier] →+* 𝒪[K.carrier]) f) (coeff_zero_map_twist _ hf0)
          (coeff_one_map_twist _ hf1) (map_residue_map_twist φ hf) := by
  rw [lubinTateTorsionSet, lubinTateTorsionSet, Set.image_iUnion]
  refine Set.iUnion_congr fun n => ?_
  rw [← Finset.coe_image,
    image_iteratedLubinTateTorsionPoints K hq φ hπmax hπne0 f hf0 hf1 hf n Φ hcompat]

def image_lubinTateLevelField.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 8, item := "Lemma 4.6", sectionId := "lemma-4-6" }

/-- ★★★★★**(i) 本体 —— `σ(K(Λ_{f,n})) = K(Λ_{f^σ,n})`**。

原文 (Yoshida08 p.8):
> Lemma 4.6. Let f, f′ ∈ O[scr]_L[X] be as above with linear coefficients π, π′, respectively. If θ ∈ Θ^L,×_π,π′ (see Corollary 3.7(ii)), then for all m ≥ 1, it gives an isomorphism [θ] = [θ]_f,f′ : µ_f,m → µ_f′,m of O[scr]-modules, and L^m_f = L^m_f′.

★★**原典と同じ主張ではない**。原典は `[θ]`(`𝒪`-線型)による輸送で
`L^m_f = L^m_{f′}` を出すが、ここでは `σ`(体の**半線型**自己同型)による輸送で
`σ(K(Λ_{f,n})) = K(Λ_{f^σ,n})` を出す。★役割(「塔が `f` の取り替えで動く」段)が
同じなので `.src` を Lemma 4.6 に付けた。逸脱として module docstring に記録した。

★★これが直前の波の申し送り (i) である。 -/
theorem image_lubinTateLevelField (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal 𝒪[K.carrier]) 𝒪[K.carrier]]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField 𝒪[K.carrier]) pp]
    [Fintype (IsLocalRing.ResidueField 𝒪[K.carrier])]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField 𝒪[K.carrier]) = pp ^ ff)
    (σ : K.carrier ≃ₐ[ℚ_[p]] K.carrier) (Φ : K.closure ≃+* K.closure)
    (hΦ : ∀ a : K.carrier, Φ (algebraMap K.carrier K.closure a)
      = algebraMap K.carrier K.closure (σ a))
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries 𝒪[K.carrier]) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) :
    Φ '' (lubinTateLevelField K hq hπmax hπne0 f hf0 hf1 hf n : Set K.closure)
      = (lubinTateLevelField K hq (maximalIdeal_eq_span_integerRingEquiv K σ hπmax)
          (integerRingEquiv_ne_zero K σ hπne0)
          (PowerSeries.map ((integerRingEquiv σ : 𝒪[K.carrier] ≃+* 𝒪[K.carrier]) :
            𝒪[K.carrier] →+* 𝒪[K.carrier]) f)
          (lubinTateSeries_integerRingEquiv K σ f hf0 hf1 hf).1
          (lubinTateSeries_integerRingEquiv K σ f hf0 hf1 hf).2.1
          (lubinTateSeries_integerRingEquiv K σ f hf0 hf1 hf).2.2 n
          : Set K.closure) := by
  rw [lubinTateLevelField, lubinTateLevelField,
    image_adjoin_of_semilinear_equiv Φ σ.toRingEquiv hΦ,
    ← Finset.coe_image,
    image_iteratedLubinTateTorsionPoints K hq (integerRingEquiv σ) hπmax hπne0 f hf0 hf1 hf n Φ
      (semilinear_integerRingEquiv K σ (Φ : K.closure →+* K.closure) hΦ)]

/-- ★★★**(i) の塔全体の版 —— `σ(K_π) = K_{σ(π)}`**(`f^σ` から作った塔)。

★**塔が一致するとは主張していない**(それが素元非依存性で、
`lubinTateClosure_sup_unramifiedClosure_eq_of_uniformizer` の仕事である)。 -/
theorem image_lubinTateClosure (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal 𝒪[K.carrier]) 𝒪[K.carrier]]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField 𝒪[K.carrier]) pp]
    [Fintype (IsLocalRing.ResidueField 𝒪[K.carrier])]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField 𝒪[K.carrier]) = pp ^ ff)
    (σ : K.carrier ≃ₐ[ℚ_[p]] K.carrier) (Φ : K.closure ≃+* K.closure)
    (hΦ : ∀ a : K.carrier, Φ (algebraMap K.carrier K.closure a)
      = algebraMap K.carrier K.closure (σ a))
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries 𝒪[K.carrier]) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f = PowerSeries.X ^ (pp ^ ff)) :
    Φ '' (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf : Set K.closure)
      = (lubinTateClosure K hq (maximalIdeal_eq_span_integerRingEquiv K σ hπmax)
          (integerRingEquiv_ne_zero K σ hπne0)
          (PowerSeries.map ((integerRingEquiv σ : 𝒪[K.carrier] ≃+* 𝒪[K.carrier]) :
            𝒪[K.carrier] →+* 𝒪[K.carrier]) f)
          (lubinTateSeries_integerRingEquiv K σ f hf0 hf1 hf).1
          (lubinTateSeries_integerRingEquiv K σ f hf0 hf1 hf).2.1
          (lubinTateSeries_integerRingEquiv K σ f hf0 hf1 hf).2.2
          : Set K.closure) := by
  rw [lubinTateClosure, lubinTateClosure,
    image_adjoin_of_semilinear_equiv Φ σ.toRingEquiv hΦ,
    image_lubinTateTorsionSet K hq (integerRingEquiv σ) hπmax hπne0 f hf0 hf1 hf Φ
      (semilinear_integerRingEquiv K σ (Φ : K.closure →+* K.closure) hΦ)]

/-! ## 7. ★★具体層 3 —— `Γ_F` の `g` を `L_S` の代数閉包へ降ろす橋

★★この節が「`F`・`S`・`g` の層」と「1 つの局所体 `K` の層」を繋ぐ。
★**#59 は出てこない**:中間体を 1 層も作らず、`closureEquivFixedField` を
`RingEquiv` に潰す(定型 (d))。 -/

/-- **★★`g` を `L_S` 自身の代数閉包の上の環自己同型として見たもの**
`Φ_g := Θ^{-1} ∘ g ∘ Θ`(`Θ := closureEquivFixedField`)。 -/
noncomputable def fixedFieldClosureAut (F : PAdicLocalField p) (S : Subgroup F.absGal)
    (hopen : IsOpen (S : Set F.absGal)) (g : F.absGal) :
    (fixedFieldLocalField F S hopen).closure ≃+* (fixedFieldLocalField F S hopen).closure :=
  ((closureEquivFixedField F S hopen).toRingEquiv.trans
      (g.toRingEquiv : F.closure ≃+* F.closure)).trans
    (closureEquivFixedField F S hopen).symm.toRingEquiv

/-- ★**`Φ_g` は `σ_g` について半線型**。★`Θ` は `L_S`-代数同型なので
`AlgEquiv.commutes` 2 回と `fixedFieldAut` の定義(`rfl`)で出る。 -/
theorem fixedFieldClosureAut_semilinear (F : PAdicLocalField p) (S : Subgroup F.absGal)
    [S.Normal] (hopen : IsOpen (S : Set F.absGal)) (g : F.absGal)
    (a : (fixedFieldLocalField F S hopen).carrier) :
    fixedFieldClosureAut F S hopen g
        (algebraMap (fixedFieldLocalField F S hopen).carrier
          (fixedFieldLocalField F S hopen).closure a)
      = algebraMap (fixedFieldLocalField F S hopen).carrier
          (fixedFieldLocalField F S hopen).closure
          (((fixedFieldAut S g).restrictScalars ℚ_[p]) a) := by
  show (closureEquivFixedField F S hopen).symm
      (g ((closureEquivFixedField F S hopen)
        (algebraMap (fixedFieldLocalField F S hopen).carrier _ a))) = _
  rw [AlgEquiv.commutes]
  have h : g (algebraMap (fixedFieldLocalField F S hopen).carrier F.closure a)
      = algebraMap (fixedFieldLocalField F S hopen).carrier F.closure
        (((fixedFieldAut S g).restrictScalars ℚ_[p]) a) := rfl
  rw [h]
  exact (closureEquivFixedField F S hopen).symm.commutes _

/-- ★`Φ_g` は整数環の上では `fixedFieldIntegerAut F S hopen g` について半線型。
★直前の波の `fixedFieldIntegerAut` の定義にちょうど当たる。 -/
theorem fixedFieldClosureAut_semilinear_integer (F : PAdicLocalField p) (S : Subgroup F.absGal)
    [S.Normal] (hopen : IsOpen (S : Set F.absGal)) (g : F.absGal)
    (a : 𝒪[(fixedFieldLocalField F S hopen).carrier]) :
    fixedFieldClosureAut F S hopen g
        (algebraMap 𝒪[(fixedFieldLocalField F S hopen).carrier]
          (fixedFieldLocalField F S hopen).closure a)
      = algebraMap 𝒪[(fixedFieldLocalField F S hopen).carrier]
          (fixedFieldLocalField F S hopen).closure (fixedFieldIntegerAut F S hopen g a) :=
  semilinear_integerRingEquiv (fixedFieldLocalField F S hopen)
    ((fixedFieldAut S g).restrictScalars ℚ_[p])
    ((fixedFieldClosureAut F S hopen g :
      (fixedFieldLocalField F S hopen).closure ≃+*
        (fixedFieldLocalField F S hopen).closure) :
      (fixedFieldLocalField F S hopen).closure →+* (fixedFieldLocalField F S hopen).closure)
    (fixedFieldClosureAut_semilinear F S hopen g) a

/-- ★★★**`Ψ_g`(`absGalConjCME`)は `Φ_g` による共役そのもの** —— `rfl` で通る。

★★これが「`Γ_F` の層」と「1 つの局所体の層」を繋ぐ橋である。
★申し送りは「#59 の (a) か (e) になる見込み」と読んでいたが、
`closureEquivFixedField` を `RingEquiv` に潰した時点で中間体の層が消え、
**定義的に等しい**ところまで落ちた(#59 の定型 (d))。 -/
theorem absGalConjCME_apply (F : PAdicLocalField p) (S : Subgroup F.absGal)
    [S.Normal] (hopen : IsOpen (S : Set F.absGal)) (g : F.absGal)
    (τ : (fixedFieldLocalField F S hopen).absGal)
    (y : (fixedFieldLocalField F S hopen).closure) :
    absGalConjCME F S hopen g τ y
      = fixedFieldClosureAut F S hopen g (τ ((fixedFieldClosureAut F S hopen g).symm y)) := rfl

/-- ★★★★**(i) を `Γ_F` の `g` に代入した形** ——
`g (L_S(Λ_{f,n})) = L_S(Λ_{f^{σ_g},n})`。★次の節点はここから持てる。 -/
theorem image_lubinTateLevelField_fixedFieldClosureAut (F : PAdicLocalField p)
    (S : Subgroup F.absGal) [S.Normal] (hopen : IsOpen (S : Set F.absGal)) (g : F.absGal)
    [IsAdicComplete (IsLocalRing.maximalIdeal 𝒪[(fixedFieldLocalField F S hopen).carrier])
      𝒪[(fixedFieldLocalField F S hopen).carrier]]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField
      𝒪[(fixedFieldLocalField F S hopen).carrier]) pp]
    [Fintype (IsLocalRing.ResidueField 𝒪[(fixedFieldLocalField F S hopen).carrier])]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField
      𝒪[(fixedFieldLocalField F S hopen).carrier]) = pp ^ ff)
    {π : 𝒪[(fixedFieldLocalField F S hopen).carrier]}
    (hπmax : IsLocalRing.maximalIdeal 𝒪[(fixedFieldLocalField F S hopen).carrier]
      = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries 𝒪[(fixedFieldLocalField F S hopen).carrier])
    (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue 𝒪[(fixedFieldLocalField F S hopen).carrier]) f
      = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) :
    fixedFieldClosureAut F S hopen g ''
        (lubinTateLevelField (fixedFieldLocalField F S hopen) hq hπmax hπne0 f hf0 hf1 hf n
          : Set (fixedFieldLocalField F S hopen).closure)
      = (lubinTateLevelField (fixedFieldLocalField F S hopen) hq
          (maximalIdeal_eq_span_integerRingEquiv (fixedFieldLocalField F S hopen)
            ((fixedFieldAut S g).restrictScalars ℚ_[p]) hπmax)
          (integerRingEquiv_ne_zero (fixedFieldLocalField F S hopen)
            ((fixedFieldAut S g).restrictScalars ℚ_[p]) hπne0)
          (PowerSeries.map ((fixedFieldIntegerAut F S hopen g :
              𝒪[(fixedFieldLocalField F S hopen).carrier] ≃+*
                𝒪[(fixedFieldLocalField F S hopen).carrier]) :
            𝒪[(fixedFieldLocalField F S hopen).carrier] →+*
              𝒪[(fixedFieldLocalField F S hopen).carrier]) f)
          (lubinTateSeries_integerRingEquiv (fixedFieldLocalField F S hopen)
            ((fixedFieldAut S g).restrictScalars ℚ_[p]) f hf0 hf1 hf).1
          (lubinTateSeries_integerRingEquiv (fixedFieldLocalField F S hopen)
            ((fixedFieldAut S g).restrictScalars ℚ_[p]) f hf0 hf1 hf).2.1
          (lubinTateSeries_integerRingEquiv (fixedFieldLocalField F S hopen)
            ((fixedFieldAut S g).restrictScalars ℚ_[p]) f hf0 hf1 hf).2.2 n
          : Set (fixedFieldLocalField F S hopen).closure) :=
  image_lubinTateLevelField (fixedFieldLocalField F S hopen) hq
    ((fixedFieldAut S g).restrictScalars ℚ_[p]) (fixedFieldClosureAut F S hopen g)
    (fixedFieldClosureAut_semilinear F S hopen g) hπmax hπne0 f hf0 hf1 hf n

end ABC3.Found.PGC
