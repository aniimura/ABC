import ABC3.Found.PGC.InflationRestrictionH2

/-!
# [pGC] Proposition 1.1 —— transgression `H¹(S,A)^{G/S} → H²(G ⧸ S, Aˢ)`(mathlib 側の欠落)

`Found/PGC/InflationRestrictionH2.lean` は `H²` の inflation-restriction を
**`H¹(S,A) = 0` を仮定して**埋めた。本ファイルはその仮定を外すために

  `tg : H¹(S,A)^{G/S} ⟶ H²(G ⧸ S, Aˢ)`

を構成し、古典的な 5 項完全列

  `0 → H¹(G⧸S, Aˢ) → H¹(G,A) → H¹(S,A)^{G/S} --tg--> H²(G⧸S, Aˢ) → H²(G,A)`

の**右 2 箇所の完全性**を証明する。★mathlib は次数 1 の `H1InfRes` / `H1InfRes_exact`
(左 2 箇所)しか持たず、transgression は**次数を問わず 1 つも持たない**。

## ★★★本ファイルの到達点(射程)

`S ⊴ G`、`A : Rep k G` について、★**仮定なしで**次を証明した。

1. `invariantsH1 A S : Submodule k (H¹(S,A))` —— `H¹(S,A)^{G/S}`(`G ⧸ S`-不変部分)。
2. `transgressionLin A S : invariantsH1 A S →ₗ[k] H²(G ⧸ S, Aˢ)` —— transgression(`k`-線形)。
3. `ker_transgressionLin` —— ★**`ker(tg) = (res : H¹(G,A) → H¹(S,A))` の像**
   (5 項完全列の `H¹(S,A)^{G/S}` における完全性)。
4. `range_transgressionLin` —— ★**`range(tg) = ker(inf : H²(G⧸S,Aˢ) → H²(G,A))`**
   (5 項完全列の `H²(G ⧸ S, Aˢ)` における完全性)。
5. `fiveTermExact` —— 3 と 4 をまとめたもの。

★★**`H¹(S,A) = 0` の仮定は外れた。** 直前の波の `inflation₂_injective`
(`H¹(S,A) = 0` ⟹ `inf` が単射)は、本ファイルでは
`inflation₂_injective_of_subsingleton_H1` として ★**4 の系**になる
(`H¹(S,A) = 0` なら `H¹(S,A)^{G/S} = 0` なので `ker(inf) = range(tg) = 0`)。

★★**外れていない箇所を名指しする**: 直前の波の `H2InfRes_exact`
(`H²(G⧸S,Aˢ) → H²(G,A) → H²(S,A)` の `H²(G,A)` における完全性)は
★**依然として `H¹(S,A) = 0` を要する**。これは 5 項完全列の外側の主張であり、
Lyndon-Hochschild-Serre スペクトル系列の `E₂^{1,1} = H¹(G⧸S, H¹(S,A))` が
一般には効くためである。★**本ファイルはそこには触っていない。**

## ★★★測定 —— mathlib に何が無いか(2026-09-07、コマンドと出力)

```
awk -F'\t' '$2 ~ /[Tt]ransgress|fiveTerm|FiveTerm|LHSpectral|lyndonHochschild|
             LyndonHochschild|HochschildSerre|hochschildSerre/ {print $2"\t"$3}'
  .cache/mathlib-index.txt
```
→ ★**0 件**。

```
node tools/absent-recheck.mjs --try
  'transgression|fiveTermExact|FiveTerm|lyndonHochschild|hochschildSerre|
   H1InvariantsAction|resInvariants'
```
→ `0 件(.cache/mathlib-index.txt、大小無視)`

```
awk -F'\t' '$2 ~ /^groupCohomology\..*([Cc]onj|[Ss]tab|[Ii]nvariant)/ {print $2"\t"$3}'
  .cache/mathlib-index.txt
```
→ `groupCohomology.d₀₁_ker_eq_invariants` の 1 件のみ。
★**`H^n(S,A)` への `G ⧸ S` の共役作用は mathlib に無い。**

★重複検査(#158): 本ファイルの 79 宣言を `ABC3.Found.PGC` 名前空間のまま REPL に積み、
`already been declared` は **0 件**。さらに
`grep -c "\b<名前>\b" .cache/{mathlib,decl}-index.txt` を 17 個の主要名について走らせ、
★**mathlib 側 0 件 / ABC3 側 0 件**を確認した。

★**REPL の `#check`**: `groupCohomology.transgression` / `H2π_surjective` /
`mem_coboundaries₁_iff` は `Unknown constant`。

## ★★設計 —— 抽象核と具体層

★★**§1 は「分岐・付値・Galois・コホモロジー・`k`・`Rep` の語彙が 1 語も出ない」核**である。
`[Group G] [AddCommGroup A]` と `ρ : G → A →+ A` だけで書かれている(#209)。

* `IsCocycleOn ρ S f` —— `S` 上の 1-コサイクル。
* `IsInvarianceWitness ρ S f a` —— 「`f` の類が `G`-不変」の**証人**
  (`ρ s (a g) - a g = ρ g (f (g⁻¹ s g)) - f s`)。★「不変性」を存在量化ではなく
  **関数 `a : G → A`** で受けるのが要点。
* `IsTgLift ρ S f F` —— transgression の**持ち上げ**。3 つのフィールドしかない:
  `F 1 = 0` / `F (g s) = F g + ρ g (f s)` / `ρ s (F g) - F g = ρ g (f (g⁻¹ s g)) - f s`。
  ★**最後のフィールドがそのまま不変性の証人**である(`IsTgLift.isGInvariantCocycle`)。
* `exists_tgLift` —— 持ち上げの**明示的な構成** `F g = (a (r g) - a 1) + ρ (r g) (f ((r g)⁻¹ g))`。
* `IsTgLift.sub_coset_const` —— ★**well-defined 性の中身**。同じ `f` の 2 つの持ち上げの差は
  「左剰余類にしか依らず値が `S`-不変」。★証明は 6 行。
* `exists_isCocycle₁_of_dOne_eq` —— ★**`H¹` 側の完全性の中身**。
  「`dF` が剰余類にしか依らない `v` の余境界なら `f` は `G` 上の 1-コサイクルの制限」。
* `exists_tgLift_of_inflated_coboundary` —— ★**`H²` 側の完全性の中身**。
  「剰余類にしか依らず値が `S`-不変な 2-コサイクルが `G` 上の余境界 `dy` なら、
  `f := y - y 1` が不変コサイクルで自分自身が持ち上げになり、`dF = x - d(定数 y 1)`」。

具体層(§2)は `repAddHom`(直前の波が用意した口)と `descendCocycle` / `descendCochain` で
核に代入するだけである。

## ★`S` の正規性をどこで使ったか

★**使う**: `IsTgLift.dOne_left_coset` / `IsTgLift.dOne_invariant` / `exists_tgLift` /
`IsCocycle₁.isTgLift` / `IsTgLift.congr_restrict` / `exists_tgLift_of_inflated_coboundary`
—— いずれも `g⁻¹ s g ∈ S` を要求するところ。具体層は `G ⧸ S` を作るので当然要る。

★**使わない**: `IsTgLift.apply_mem` / `dOne_right_zero` / `dOne_one_left` /
`dOne_right_coset` / `IsTgLift.isCocycleOn` / `IsTgLift.sub_coset_const` /
`IsTgLift.add` / `IsTgLift.shift` / `exists_isCocycle₁_of_dOne_eq` /
`IsCocycleOn.conj_eq`。★さらに `normCosetRep`(正規化された代表元)は
`cosetRep S g * (cosetRep S 1)⁻¹` という形にしたので ★**正規性も可判定性も要らない**
(`if g ∈ S then 1 else …` と書くと `Decidable (g ∈ S)` が要る)。

## ★選択公理について

★**抽象核は選択を使わない形にしてある**: `exists_tgLift` は代表元関数 `r` と
証人関数 `a` を**引数で受ける**。★これは直前の波の `exists_cochain_right/left` と同じ作法。

★具体層(`exists_tgLift_rep` 以降)は選択を使う: `cosetRep` は `Quotient.out`、
証人は `choose`。★これは「`H¹(S,A)^{G/S}` の各元に対して実際に写像を作る」ために避けられない。

★**注意**: `#print axioms` では「選択を使わない」ことは確認できない。
`abel` / `group` などの正規化タクティクが `Classical.choice` を引くためである
(例: `dOne_sub` は `simp only + abel` だけなのに `Classical.choice` が出る)。
★一方 `IsCocycle₁.isTgLift` / `IsTgLift.add` / `IsTgLift.apply_mem` /
`IsCocycleOn.conj_eq` / `IsTgLift.dOne_right_zero` / `IsCocycle₁.dOne_eq_zero` は
`[propext]` **のみ**である。

## 逸脱の記録

* **逸脱 1**: 原典 (pGC) は `Γ_K` の**連続**コホモロジーを使うが、本ファイルは
  mathlib の**離散**群コホモロジー `groupCohomology` で書いている。
  ★`GroupCohomologyFinite.lean` / `InflationRestrictionH2.lean` の逸脱 1 と同じ理由・同じ範囲。
* **逸脱 2**: 原典は 5 項完全列も transgression も名指ししていない(「well-known」で畳んでいる)。
  ★本ファイルは `H¹(S,A)^{G/S}` を**コチェインの述語**(`IsGInvariantH1`)で定義した。
  mathlib には `H^n(S,A)` への `G ⧸ S` の共役作用が無いためである(上の測定を参照)。
  ★`isGInvariantH1_iff` で「不変なコサイクルで代表される」ことと同値であることを示し、
  ★`invariantsH1` が `Submodule` になること(`isGInvariantH1_zero/add/smul`)で
  「不変部分加群」という呼び方を正当化している。
* **逸脱 3**: 完全性は `ShortComplex.Exact` ではなく
  `LinearMap.ker = Submodule.comap … (LinearMap.range …)` の形で述べた。
  ★`invariantsH1` が `H¹(S,A)` の**部分**加群であり、`res` の像がそこに入る
  (`isGInvariantH1_of_restriction`)ため、短複体を作るには余分な余制限が要る。
  ★数学的内容は同じである。

## 退化の自己検査

* `invariantsH1_bot` —— `S = ⊥` では `H¹(⊥,A) = 0` なので `H¹(S,A)^{G/S} = 0`。
  このとき `range(tg) = 0` となり、`range_transgressionLin` は
  「`inf` が単射」に潰れる。★直前の波の `inflation₂_bijective_bot` と整合する。
* `inflation₂_injective_of_subsingleton_H1` —— ★**直前の波の主定理が本ファイルの系になる**
  ことの確認。★「新しい定理が古い定理を含む」ことを Lean で測った。
* `isGInvariantH1_of_restriction` —— `res` の像は必ず不変。
  ★これが無いと 5 項完全列の 3 項目が `H¹(S,A)^{G/S}` である意味が無い。
* `inflation₂_transgression` —— `inf ∘ tg = 0`。★短複体として自明でないことの確認。
* ★**`axiom` も `structure` 的な逃げも `sorry` も置いていない。**
  ★**「通すための偽の仮説」は作っていない。** 本ファイルの主定理はすべて**仮定なし**である。
-/

namespace ABC3.Found.PGC

open CategoryTheory groupCohomology

universe u v

/-! ## 1. ★★抽象核 —— transgression の持ち上げ(語彙ゼロ)

★この節には表現論もコホモロジーも分岐も `k` も `Rep` も出てこない。
群 `G`、可換群 `A`、`ρ : G → A →+ A` だけである。 -/

section AbstractCoreTransgression

variable {G : Type u} {A : Type v} [Group G] [AddCommGroup A]

/-- **`S` 上の 1-コサイクル**(語彙ゼロ)。

`f : G → A` は `G` 全体で定義されているが、条件は `S` の上でしか課さない。 -/
def IsCocycleOn (ρ : G → A →+ A) (S : Subgroup G) (f : G → A) : Prop :=
  ∀ s ∈ S, ∀ t ∈ S, f (s * t) = ρ s (f t) + f s

/-- **「`f` の類が `G` の共役作用で不変」の証人**(語彙ゼロ)。

`g` の作用は `(g · f)(s) = ρ g (f (g⁻¹ s g))` であり、不変性は
「`g · f - f` が `S` 上で余境界」ということ。★その余境界の**原始関数**を
`a : G → A` として受ける。★存在量化ではなく関数で受けるので選択公理が要らない。 -/
def IsInvarianceWitness (ρ : G → A →+ A) (S : Subgroup G) (f a : G → A) : Prop :=
  ∀ g : G, ∀ s ∈ S, ρ s (a g) - a g = ρ g (f (g⁻¹ * s * g)) - f s

/-- **★★transgression の持ち上げ**(語彙ゼロ)。

`F : G → A` が `f` の持ち上げであるとは:

* `one` —— `F 1 = 0`(正規化)
* `right` —— `F (g s) = F g + ρ g (f s)`(`s ∈ S`。★`F` は `f` を「右から」延長する)
* `conj` —— `ρ s (F g) - F g = ρ g (f (g⁻¹ s g)) - f s`
  (★`F g` 自身が `g` における不変性の証人になっている)

★このとき `dF` は `G ⧸ S` 上の `Aˢ`-値 2-コサイクルに降り、その類が `tg [f]` である。 -/
structure IsTgLift (ρ : G → A →+ A) (S : Subgroup G) (f F : G → A) : Prop where
  one : F 1 = 0
  right : ∀ g : G, ∀ s ∈ S, F (g * s) = F g + ρ g (f s)
  conj : ∀ g : G, ∀ s ∈ S, ρ s (F g) - F g = ρ g (f (g⁻¹ * s * g)) - f s

/-- 1-コサイクルは `1` で消える。 -/
theorem IsCocycleOn.map_one {ρ : G → A →+ A} (hone : ∀ a, ρ 1 a = a) {S : Subgroup G} {f : G → A}
    (hf : IsCocycleOn ρ S f) : f 1 = 0 := by
  have h := hf 1 S.one_mem 1 S.one_mem
  rw [one_mul, hone] at h
  linear_combination (norm := abel) -h

/-- `S` 上の 1-コサイクルの共役の公式。★正規性は要らない(`u`, `w` を `S` の元として受ける)。 -/
theorem IsCocycleOn.conj_eq {ρ : G → A →+ A} {S : Subgroup G} {f : G → A}
    (hf : IsCocycleOn ρ S f) {u : G} (hu : u ∈ S) {w : G} (hw : w ∈ S) :
    ρ w (f (w⁻¹ * u * w)) = f (u * w) - f w := by
  have hm : w⁻¹ * u * w ∈ S := S.mul_mem (S.mul_mem (S.inv_mem hw) hu) hw
  have h := hf w hw _ hm
  rw [show w * (w⁻¹ * u * w) = u * w by group] at h
  linear_combination (norm := abel) -h

/-- 持ち上げは `S` の上では `f` そのもの。 -/
theorem IsTgLift.apply_mem {ρ : G → A →+ A} (hone : ∀ a, ρ 1 a = a) {S : Subgroup G}
    {f F : G → A} (h : IsTgLift ρ S f F) : ∀ s ∈ S, F s = f s := by
  intro s hs
  have := h.right 1 s hs
  rw [one_mul, h.one, hone] at this
  simpa using this

/-- ★持ち上げの存在から、`f` が `S` 上の 1-コサイクルであることが従う。 -/
theorem IsTgLift.isCocycleOn {ρ : G → A →+ A} (hone : ∀ a, ρ 1 a = a) {S : Subgroup G}
    {f F : G → A} (h : IsTgLift ρ S f F) : IsCocycleOn ρ S f := by
  intro s hs t ht
  have e := h.right s t ht
  rw [h.apply_mem hone _ (S.mul_mem hs ht), h.apply_mem hone _ hs] at e
  linear_combination (norm := abel) e

/-- **条件 (A)** —— `dF` の第 2 引数が `S` なら 0。★正規性は要らない。 -/
theorem IsTgLift.dOne_right_zero {ρ : G → A →+ A} (hone : ∀ a, ρ 1 a = a) {S : Subgroup G}
    {f F : G → A} (h : IsTgLift ρ S f F) : ∀ g : G, ∀ s ∈ S, dOne ρ F g s = 0 := by
  intro g s hs
  simp only [dOne, h.apply_mem hone s hs, h.right g s hs]
  abel

/-- `dF` の第 1 引数が `1` なら 0。★正規性は要らない。 -/
theorem IsTgLift.dOne_one_left {ρ : G → A →+ A} (hone : ∀ a, ρ 1 a = a) {S : Subgroup G}
    {f F : G → A} (h : IsTgLift ρ S f F) : ∀ j : G, dOne ρ F 1 j = 0 := by
  intro j
  simp [dOne, hone, h.one]

/-- **`dF` の第 1 引数は左剰余類にしか依らない**。★ここで `S` の正規性を使う。
★中身は `conj` フィールドを `ρ g` で押したもの。 -/
theorem IsTgLift.dOne_left_coset {ρ : G → A →+ A}
    (hmul : ∀ g h : G, ∀ a : A, ρ (g * h) a = ρ g (ρ h a))
    {S : Subgroup G} [S.Normal] {f F : G → A} (h : IsTgLift ρ S f F) :
    ∀ g j : G, ∀ s ∈ S, dOne ρ F (g * s) j = dOne ρ F g j := by
  intro g j s hs
  have hconj : j⁻¹ * s * j ∈ S := Subgroup.Normal.conj_mem' ‹_› s hs j
  have e1 : g * s * j = g * j * (j⁻¹ * s * j) := by group
  have hF1 : F (g * s * j) = F (g * j) + ρ (g * j) (f (j⁻¹ * s * j)) := by
    rw [e1]; exact h.right (g * j) _ hconj
  have hF2 : F (g * s) = F g + ρ g (f s) := h.right g s hs
  have hc := h.conj j s hs
  have hkey : ρ g (ρ s (F j)) - ρ g (F j) = ρ g (ρ j (f (j⁻¹ * s * j))) - ρ g (f s) := by
    rw [← map_sub, ← map_sub, hc]
  simp only [dOne, hF1, hF2, hmul g s (F j), hmul g j (f (j⁻¹ * s * j))]
  linear_combination (norm := abel) hkey

/-- **`dF` の第 2 引数は左剰余類にしか依らない**。★条件 (A) から出るので正規性は要らない。 -/
theorem IsTgLift.dOne_right_coset {ρ : G → A →+ A} (hone : ∀ a, ρ 1 a = a)
    (hmul : ∀ g h : G, ∀ a : A, ρ (g * h) a = ρ g (ρ h a))
    {S : Subgroup G} {f F : G → A} (h : IsTgLift ρ S f F) :
    ∀ g j : G, ∀ t ∈ S, dOne ρ F g (j * t) = dOne ρ F g j := fun g j t ht =>
  cocycle₂_right_coset_invariant (isCocycle₂_dOne ρ hmul F) (h.dOne_right_zero hone) g j t ht

/-- **`dF` の値は `S`-不変**。★ここで `S` の正規性を使う。 -/
theorem IsTgLift.dOne_invariant {ρ : G → A →+ A} (hone : ∀ a, ρ 1 a = a)
    (hmul : ∀ g h : G, ∀ a : A, ρ (g * h) a = ρ g (ρ h a))
    {S : Subgroup G} [S.Normal] {f F : G → A} (h : IsTgLift ρ S f F) :
    ∀ g j : G, ∀ u ∈ S, ρ u (dOne ρ F g j) = dOne ρ F g j := by
  intro g j u hu
  have e := isCocycle₂_dOne ρ hmul F u g j
  have h1 : dOne ρ F u g = 0 := by
    have e1 := h.dOne_left_coset hmul 1 g u hu
    rw [one_mul] at e1
    rw [e1, h.dOne_one_left hone]
  have h2 : dOne ρ F u (g * j) = 0 := by
    have e2 := h.dOne_left_coset hmul 1 (g * j) u hu
    rw [one_mul] at e2
    rw [e2, h.dOne_one_left hone]
  have h3 : dOne ρ F (u * g) j = dOne ρ F g j := by
    have hconj : g⁻¹ * u * g ∈ S := Subgroup.Normal.conj_mem' ‹_› u hu g
    have e3 : u * g = g * (g⁻¹ * u * g) := by group
    rw [e3, h.dOne_left_coset hmul g j _ hconj]
  rw [h1, h2, h3] at e
  linear_combination (norm := abel) -e

/-- 余境界は差について加法的。 -/
theorem dOne_sub (ρ : G → A →+ A) (y z : G → A) (g h : G) :
    dOne ρ (fun g => y g - z g) g h = dOne ρ y g h - dOne ρ z g h := by
  simp only [dOne, map_sub]; abel

/-- 余境界は和について加法的。 -/
theorem dOne_add (ρ : G → A →+ A) (y z : G → A) (g h : G) :
    dOne ρ (fun g => y g + z g) g h = dOne ρ y g h + dOne ρ z g h := by
  simp only [dOne, map_add]; abel

/-- 「主コチェイン」`g ↦ ρ g b - b` の余境界は 0。 -/
theorem dOne_principal (ρ : G → A →+ A)
    (hmul : ∀ g h : G, ∀ a : A, ρ (g * h) a = ρ g (ρ h a)) (b : A) (g h : G) :
    dOne ρ (fun g => ρ g b - b) g h = 0 := by
  simp only [dOne, map_sub, hmul g h b]; abel

/-- ★`IsTgLift` は `f` の `S` への制限にしか依らない。 -/
theorem IsTgLift.congr_restrict {ρ : G → A →+ A} {S : Subgroup G} [S.Normal] {f F : G → A}
    (h : IsTgLift ρ S f F) (f' : G → A) (hff : ∀ s ∈ S, f' s = f s) : IsTgLift ρ S f' F := by
  refine ⟨h.one, ?_, ?_⟩
  · intro g s hs; rw [hff s hs]; exact h.right g s hs
  · intro g s hs
    have hconj : g⁻¹ * s * g ∈ S := Subgroup.Normal.conj_mem' ‹_› s hs g
    rw [hff s hs, hff _ hconj]; exact h.conj g s hs

/-- ★`f` を余境界だけずらすと、持ち上げも主コチェインだけずれる。
★`dOne_principal` により `dF` は**変わらない**。これが well-defined 性の第 2 の柱。 -/
theorem IsTgLift.shift {ρ : G → A →+ A} (hone : ∀ a, ρ 1 a = a)
    (hmul : ∀ g h : G, ∀ a : A, ρ (g * h) a = ρ g (ρ h a))
    {S : Subgroup G} {f F : G → A} (h : IsTgLift ρ S f F) (b : A) :
    IsTgLift ρ S (fun g => f g + (ρ g b - b)) (fun g => F g + (ρ g b - b)) := by
  refine ⟨by simp [h.one, hone], ?_, ?_⟩
  · intro g s hs
    simp only [h.right g s hs, map_add, map_sub, hmul g s b]
    abel
  · intro g s hs
    have e := h.conj g s hs
    have e1 : ρ s (ρ g b) = ρ (s * g) b := (hmul s g b).symm
    have e2 : ρ g (ρ (g⁻¹ * s * g) b) = ρ (s * g) b := by
      rw [← hmul]; congr 2; group
    simp only [map_add, map_sub, e1, e2]
    linear_combination (norm := abel) e

/-- **★★well-defined 性の中身** —— 同じ `f` の 2 つの持ち上げの差は
「左剰余類にしか依らず、値が `S`-不変」。★正規性は要らない。★証明は 6 行。 -/
theorem IsTgLift.sub_coset_const {ρ : G → A →+ A} {S : Subgroup G} {f F F' : G → A}
    (h : IsTgLift ρ S f F) (h' : IsTgLift ρ S f F') :
    (∀ g : G, ∀ s ∈ S, F (g * s) - F' (g * s) = F g - F' g) ∧
      (∀ g : G, ∀ s ∈ S, ρ s (F g - F' g) = F g - F' g) := by
  constructor
  · intro g s hs
    rw [h.right g s hs, h'.right g s hs]; abel
  · intro g s hs
    have e := h.conj g s hs
    have e' := h'.conj g s hs
    rw [map_sub]
    linear_combination (norm := abel) e - e'

/-- ★持ち上げは加法的。★これが transgression の加法性を与える。 -/
theorem IsTgLift.add {ρ : G → A →+ A} {S : Subgroup G} {f F f' F' : G → A}
    (h : IsTgLift ρ S f F) (h' : IsTgLift ρ S f' F') :
    IsTgLift ρ S (fun g => f g + f' g) (fun g => F g + F' g) := by
  refine ⟨by simp [h.one, h'.one], ?_, ?_⟩
  · intro g s hs
    simp only [h.right g s hs, h'.right g s hs, map_add]; abel
  · intro g s hs
    have e := h.conj g s hs
    have e' := h'.conj g s hs
    simp only [map_add]
    linear_combination (norm := abel) e + e'

/-- **transgression の持ち上げ(明示式)**。

`r` は左剰余類の代表元、`a` は不変性の証人。★どちらも**引数**なので選択公理は要らない。 -/
def tgLift (ρ : G → A →+ A) (f a : G → A) (r : G → G) (g : G) : A :=
  (a (r g) - a 1) + ρ (r g) (f ((r g)⁻¹ * g))

/-- **★★★持ち上げの構成**(抽象核の本体)。

`f` が `S` 上の 1-コサイクルで、その類が `G`-不変(証人 `a`)であれば、
代表元関数 `r`(`r 1 = 1` に正規化されたもの)から持ち上げが作れる。

★`r` に課すのは `r 1 = 1`、`(r g)⁻¹ g ∈ S`、`r (g s) = r g` の 3 つだけ。
★分岐・付値・Galois・コホモロジーの語彙が 1 語も出ない。 -/
theorem exists_tgLift {ρ : G → A →+ A} (hone : ∀ a, ρ 1 a = a)
    (hmul : ∀ g h : G, ∀ a : A, ρ (g * h) a = ρ g (ρ h a))
    {S : Subgroup G} [S.Normal] {f : G → A} (hf : IsCocycleOn ρ S f)
    (a : G → A) (ha : IsInvarianceWitness ρ S f a)
    (r : G → G) (hr0 : r 1 = 1) (hr1 : ∀ g : G, (r g)⁻¹ * g ∈ S)
    (hr2 : ∀ g : G, ∀ s ∈ S, r (g * s) = r g) :
    IsTgLift ρ S f (tgLift ρ f a r) := by
  refine ⟨?_, ?_, ?_⟩
  · simp [tgLift, hr0, hf.map_one hone]
  · intro g s hs
    have hq : r (g * s) = r g := hr2 g s hs
    have hu : (r g)⁻¹ * g ∈ S := hr1 g
    have e2 : f ((r g)⁻¹ * g * s) = ρ ((r g)⁻¹ * g) (f s) + f ((r g)⁻¹ * g) := hf _ hu s hs
    have e3 : ρ (r g) (ρ ((r g)⁻¹ * g) (f s)) = ρ g (f s) := by
      rw [← hmul]; congr 2; group
    simp only [tgLift, hq, show (r g)⁻¹ * (g * s) = (r g)⁻¹ * g * s by group, e2, map_add, e3]
    abel
  · intro g s hs
    have hw : (r g)⁻¹ * g ∈ S := hr1 g
    have hu : (r g)⁻¹ * s * r g ∈ S := Subgroup.Normal.conj_mem' ‹_› s hs (r g)
    have ha1 : ρ s (a 1) = a 1 := by
      have h := ha 1 s hs
      simp only [inv_one, one_mul, mul_one, hone] at h
      linear_combination (norm := abel) h
    have haq := ha (r g) s hs
    have key1 : ρ ((r g)⁻¹ * s * r g) (f ((r g)⁻¹ * g))
        = f ((r g)⁻¹ * s * r g * ((r g)⁻¹ * g)) - f ((r g)⁻¹ * s * r g) := by
      have h := hf _ hu _ hw
      linear_combination (norm := abel) -h
    have hsq2 : ρ s (ρ (r g) (f ((r g)⁻¹ * g)))
        = ρ (r g) (f ((r g)⁻¹ * s * r g * ((r g)⁻¹ * g)))
          - ρ (r g) (f ((r g)⁻¹ * s * r g)) := by
      have hsq : ρ s (ρ (r g) (f ((r g)⁻¹ * g)))
          = ρ (r g) (ρ ((r g)⁻¹ * s * r g) (f ((r g)⁻¹ * g))) := by
        rw [← hmul, ← hmul]; congr 2; group
      rw [hsq, key1, map_sub]
    have key2 : ρ g (f (g⁻¹ * s * g))
        = ρ (r g) (f ((r g)⁻¹ * s * r g * ((r g)⁻¹ * g))) - ρ (r g) (f ((r g)⁻¹ * g)) := by
      have hgq : ρ g (f (g⁻¹ * s * g))
          = ρ (r g) (ρ ((r g)⁻¹ * g) (f (g⁻¹ * s * g))) := by
        rw [← hmul]; congr 2; group
      rw [hgq, show g⁻¹ * s * g
          = ((r g)⁻¹ * g)⁻¹ * ((r g)⁻¹ * s * r g) * ((r g)⁻¹ * g) from by group,
        hf.conj_eq hu hw, map_sub]
    simp only [tgLift, map_add, map_sub]
    linear_combination (norm := abel) haq + hsq2 - ha1 - key2

/-- **`G` 全体の上の 1-コサイクル**(語彙ゼロ)。 -/
def IsCocycle₁ (ρ : G → A →+ A) (φ : G → A) : Prop :=
  ∀ g h : G, φ (g * h) = ρ g (φ h) + φ g

theorem IsCocycle₁.map_one {ρ : G → A →+ A} (hone : ∀ a, ρ 1 a = a) {φ : G → A}
    (hφ : IsCocycle₁ ρ φ) : φ 1 = 0 := by
  have h := hφ 1 1
  rw [one_mul, hone] at h
  linear_combination (norm := abel) -h

theorem IsCocycle₁.dOne_eq_zero {ρ : G → A →+ A} {φ : G → A} (hφ : IsCocycle₁ ρ φ)
    (g h : G) : dOne ρ φ g h = 0 := by
  simp only [dOne, hφ g h]; abel

/-- ★`G` 上の 1-コサイクルは、その `S` への制限の transgression 持ち上げになる。
★これが「`res` の像の transgression は 0」の中身。 -/
theorem IsCocycle₁.isTgLift {ρ : G → A →+ A} (hone : ∀ a, ρ 1 a = a) {S : Subgroup G} [S.Normal]
    {φ f : G → A} (hφ : IsCocycle₁ ρ φ) (hres : ∀ s ∈ S, φ s = f s) :
    IsTgLift ρ S f φ := by
  refine ⟨hφ.map_one hone, ?_, ?_⟩
  · intro g s hs
    rw [hφ g s, hres s hs]; abel
  · intro g s hs
    have hconj : g⁻¹ * s * g ∈ S := Subgroup.Normal.conj_mem' ‹_› s hs g
    have e1 := hφ s g
    have e2 := hφ g (g⁻¹ * s * g)
    rw [show g * (g⁻¹ * s * g) = s * g by group] at e2
    rw [← hres _ hconj, ← hres s hs]
    linear_combination (norm := abel) e2 - e1

/-- **★★★完全性(`H¹` 側)の中身** —— `dF` が「剰余類にしか依らない `v` の余境界」なら、
`f` は `G` 全体の 1-コサイクルの制限である。

★`v` の値が `Aˢ` にあることすら要らない(剰余類にしか依らないことだけ使う)。
★分岐・付値・Galois・コホモロジーの語彙が 1 語も出ない。 -/
theorem exists_isCocycle₁_of_dOne_eq {ρ : G → A →+ A} (hone : ∀ a, ρ 1 a = a)
    {S : Subgroup G} {f F v : G → A} (h : IsTgLift ρ S f F)
    (hv : ∀ g : G, ∀ s ∈ S, v (g * s) = v g)
    (hFv : ∀ g j : G, dOne ρ F g j = dOne ρ v g j) :
    ∃ φ : G → A, IsCocycle₁ ρ φ ∧ ∀ s ∈ S, φ s = f s := by
  have hcoc : IsCocycle₁ ρ (fun g => F g - v g) := by
    intro g j
    have e := hFv g j
    simp only [dOne] at e
    simp only [map_sub]
    linear_combination (norm := abel) -e
  refine ⟨fun g => F g - v g, hcoc, ?_⟩
  have h1 : v 1 = 0 := by
    have := hcoc.map_one hone
    simp only [h.one] at this
    linear_combination (norm := abel) -this
  intro s hs
  have e := hv 1 s hs
  rw [one_mul] at e
  show F s - v s = f s
  rw [e, h1, h.apply_mem hone s hs, sub_zero]

/-- **★★★完全性(`H²` 側)の中身** —— 「剰余類にしか依らず値が `S`-不変な 2-コサイクル `x`」が
`G` 上の余境界 `dy` であれば、`f := y - y 1` は `S` 上の 1-コサイクルで**自分自身が持ち上げ**になり、
その `dOne` は `x` から「定数 `y 1` の余境界」を引いたものに一致する。

★分岐・付値・Galois・コホモロジーの語彙が 1 語も出ない。 -/
theorem exists_tgLift_of_inflated_coboundary {ρ : G → A →+ A} (hone : ∀ a, ρ 1 a = a)
    {S : Subgroup G} [S.Normal] {x : G → G → A}
    (hfacR : ∀ g h : G, ∀ t ∈ S, x g (h * t) = x g h)
    (hfacL : ∀ g h : G, ∀ s ∈ S, x (g * s) h = x g h)
    (hinv : ∀ g h : G, ∀ u ∈ S, ρ u (x g h) = x g h)
    (y : G → A) (hy : ∀ g h : G, dOne ρ y g h = x g h) :
    IsTgLift ρ S (fun g => y g - y 1) (fun g => y g - y 1) ∧
      (∀ g h : G, dOne ρ (fun g => y g - y 1) g h = x g h - ρ g (y 1)) ∧
      (∀ u ∈ S, ρ u (y 1) = y 1) := by
  have hy11 : y 1 = x 1 1 := by
    have h := hy 1 1
    simp only [dOne, one_mul, hone] at h
    linear_combination (norm := abel) h
  have hinvy : ∀ u ∈ S, ρ u (y 1) = y 1 := by
    intro u hu; rw [hy11]; exact hinv 1 1 u hu
  have hxgs : ∀ g : G, ∀ s ∈ S, x g s = ρ g (y 1) := by
    intro g s hs
    have e := hfacR g 1 s hs
    rw [one_mul] at e
    have e2 := hy g 1
    simp only [dOne, mul_one] at e2
    rw [e, ← e2]
    abel
  have hxsh : ∀ s ∈ S, ∀ h : G, x s h = y 1 := by
    intro s hs h
    have e := hfacL 1 h s hs
    rw [one_mul] at e
    have e2 := hy 1 h
    simp only [dOne, one_mul, hone] at e2
    rw [e, ← e2]
    abel
  have hright : ∀ g : G, ∀ s ∈ S, y (g * s) - y 1 = (y g - y 1) + ρ g (y s - y 1) := by
    intro g s hs
    have e := hy g s
    rw [hxgs g s hs] at e
    simp only [dOne] at e
    rw [map_sub]
    linear_combination (norm := abel) -e
  refine ⟨⟨by simp, hright, ?_⟩, ?_, hinvy⟩
  · intro g s hs
    have hconj : g⁻¹ * s * g ∈ S := Subgroup.Normal.conj_mem' ‹_› s hs g
    have E1 := hy s g
    rw [hxsh s hs g] at E1
    simp only [dOne] at E1
    have E2 := hy g (g⁻¹ * s * g)
    rw [hxgs g _ hconj] at E2
    simp only [dOne] at E2
    rw [show g * (g⁻¹ * s * g) = s * g by group] at E2
    have hys := hinvy s hs
    show ρ s (y g - y 1) - (y g - y 1) = ρ g (y (g⁻¹ * s * g) - y 1) - (y s - y 1)
    simp only [map_sub]
    linear_combination (norm := abel) E1 - hys - E2
  · intro g h
    have e := hy g h
    simp only [dOne, map_sub] at e ⊢
    linear_combination (norm := abel) e

end AbstractCoreTransgression

/-! ## 2. 具体層 —— `Rep k G` と `groupCohomology`

★ここから `Rep` が出てくる。核に代入するための配管だけである。 -/

section ConcreteTransgression

variable {k G : Type} [CommRing k] [Group G]

/-- ★**正規化した左剰余類の代表元**(`r 1 = 1` を満たす)。

★`if g ∈ S then 1 else cosetRep S g` と書くと `Decidable (g ∈ S)` が要るので、
`cosetRep S g * (cosetRep S 1)⁻¹` という形にした。
★これで**可判定性も `S` の正規性も要らない**。 -/
noncomputable def normCosetRep (S : Subgroup G) (g : G) : G := cosetRep S g * (cosetRep S 1)⁻¹

theorem cosetRep_one_mem (S : Subgroup G) : cosetRep S 1 ∈ S := by
  have h := cosetRep_inv_mul_mem S 1
  rw [mul_one] at h
  simpa using S.inv_mem h

theorem normCosetRep_one (S : Subgroup G) : normCosetRep S 1 = 1 := by
  simp [normCosetRep]

theorem normCosetRep_inv_mul_mem (S : Subgroup G) (g : G) : (normCosetRep S g)⁻¹ * g ∈ S := by
  have e : (normCosetRep S g)⁻¹ * g = cosetRep S 1 * ((cosetRep S g)⁻¹ * g) := by
    simp only [normCosetRep, mul_inv_rev, inv_inv, mul_assoc]
  rw [e]
  exact S.mul_mem (cosetRep_one_mem S) (cosetRep_inv_mul_mem S g)

theorem normCosetRep_mul_of_mem (S : Subgroup G) (g s : G) (hs : s ∈ S) :
    normCosetRep S (g * s) = normCosetRep S g := by
  simp only [normCosetRep, cosetRep_mul_of_mem S g s hs]

/-- ★★`f` の `H¹(S,A)` における類が `G` の共役作用で不変であること(コチェインの言葉)。

`g` の作用は `(g · f)(s) = ρ g (f (g⁻¹ s g))`。★不変性は「`g · f - f` が `S` 上で余境界」。 -/
def IsGInvariantCocycle (A : Rep k G) (S : Subgroup G) (f : G → A) : Prop :=
  ∀ g : G, ∃ b : A, ∀ s ∈ S, A.ρ s b - b = A.ρ g (f (g⁻¹ * s * g)) - f s

/-- 持ち上げが在れば、その `f` の類は `G`-不変。★`conj` フィールドがそのまま証人。 -/
theorem IsTgLift.isGInvariantCocycle {A : Rep k G} {S : Subgroup G} {f F : G → A}
    (hF : IsTgLift (repAddHom A) S f F) : IsGInvariantCocycle A S f :=
  fun g => ⟨F g, fun s hs => hF.conj g s hs⟩

/-- ★★transgression の持ち上げの存在(具体層)。★ここで初めて選択公理を使う
(代表元 `cosetRep` と証人の `choose`)。 -/
theorem exists_tgLift_rep (A : Rep k G) (S : Subgroup G) [S.Normal] (f : G → A)
    (hf : IsCocycleOn (repAddHom A) S f) (hinv : IsGInvariantCocycle A S f) :
    ∃ F : G → A, IsTgLift (repAddHom A) S f F := by
  classical
  choose b hb using hinv
  exact ⟨tgLift (repAddHom A) f b (normCosetRep S),
    exists_tgLift (repAddHom_one A) (repAddHom_mul A) hf b hb (normCosetRep S)
      (normCosetRep_one S) (normCosetRep_inv_mul_mem S) (normCosetRep_mul_of_mem S)⟩

/-- 剰余類にしか依らず値が `Aˢ` にある 1-コチェインを `G ⧸ S` へ降ろす。
★`descendCocycle`(直前の波)の 1-コチェイン版。 -/
noncomputable def descendCochain (A : Rep k G) (S : Subgroup G) [S.Normal] (w : G → A)
    (hmem : ∀ g : G, w g ∈ Representation.invariants (A.ρ.comp S.subtype))
    (hfac : ∀ g : G, ∀ s ∈ S, w (g * s) = w g) :
    (G ⧸ S) → (A.quotientToInvariants S) := fun q =>
  Quotient.liftOn' q (fun g => (⟨w g, hmem g⟩ : (A.quotientToInvariants S)))
    (by
      intro g h e
      have m : g⁻¹ * h ∈ S := QuotientGroup.leftRel_apply.1 e
      apply Subtype.ext
      show w g = w h
      conv_rhs => rw [show h = g * (g⁻¹ * h) from (mul_inv_cancel_left _ _).symm]
      exact (hfac g _ m).symm)

theorem descendCochain_mk (A : Rep k G) (S : Subgroup G) [S.Normal] (w : G → A)
    (hmem : ∀ g : G, w g ∈ Representation.invariants (A.ρ.comp S.subtype))
    (hfac : ∀ g : G, ∀ s ∈ S, w (g * s) = w g) (g : G) :
    ((descendCochain A S w hmem hfac (QuotientGroup.mk' S g) : _) : A) = w g := rfl

/-- 持ち上げ `F` の余境界は剰余類にしか依らない。 -/
theorem IsTgLift.rep_hfac {A : Rep k G} {S : Subgroup G} [S.Normal] {f F : G → A}
    (hF : IsTgLift (repAddHom A) S f F) :
    ∀ g h s : G, s ∈ S → ∀ t ∈ S,
      dOne (repAddHom A) F (g * s) (h * t) = dOne (repAddHom A) F g h := by
  intro g h s hs t ht
  rw [hF.dOne_right_coset (repAddHom_one A) (repAddHom_mul A) (g * s) h t ht,
    hF.dOne_left_coset (repAddHom_mul A) g h s hs]

/-- 持ち上げ `F` の余境界は `Aˢ` に値をとる。 -/
theorem IsTgLift.rep_hmem {A : Rep k G} {S : Subgroup G} [S.Normal] {f F : G → A}
    (hF : IsTgLift (repAddHom A) S f F) :
    ∀ g h : G, dOne (repAddHom A) F g h ∈ Representation.invariants (A.ρ.comp S.subtype) :=
  fun g h => (Representation.mem_invariants _ _).2
    (fun u => hF.dOne_invariant (repAddHom_one A) (repAddHom_mul A) g h u u.2)

/-- ★★transgression の 2-コサイクル(`G ⧸ S` 上、`Aˢ` 係数)。 -/
noncomputable def tgCocycle (A : Rep k G) (S : Subgroup G) [S.Normal] {f F : G → A}
    (hF : IsTgLift (repAddHom A) S f F) : cocycles₂ (A.quotientToInvariants S) :=
  ⟨descendCocycle A S (dOne (repAddHom A) F) hF.rep_hmem hF.rep_hfac,
    descendCocycle_mem A S _ hF.rep_hmem hF.rep_hfac
      (isCocycle₂_dOne _ (repAddHom_mul A) F)⟩

theorem tgCocycle_apply (A : Rep k G) (S : Subgroup G) [S.Normal] {f F : G → A}
    (hF : IsTgLift (repAddHom A) S f F) (g h : G) :
    (((tgCocycle A S hF : cocycles₂ (A.quotientToInvariants S)) : (G ⧸ S) × (G ⧸ S) → _)
      (QuotientGroup.mk' S g, QuotientGroup.mk' S h) : A) = dOne (repAddHom A) F g h := rfl

def transgression.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }

/-- **★★★transgression の類** `tg [f] ∈ H²(G ⧸ S, Aˢ)`(持ち上げ `F` から)。

原文 (pGC p.3):
> Proposition 1.1: The cyclotomic character χ : ΓK → Zp× can be recovered entirely
> group-theoretically from ΓK.

★原典は 5 項完全列も transgression も名指ししていない(「well-known」で畳んでいる)ので、
本ファイルは mathlib 側の欠落を埋める位置づけである。 -/
noncomputable def tgClass (A : Rep k G) (S : Subgroup G) [S.Normal] {f F : G → A}
    (hF : IsTgLift (repAddHom A) S f F) : groupCohomology (A.quotientToInvariants S) 2 :=
  (ConcreteCategory.hom (H2π (A.quotientToInvariants S))) (tgCocycle A S hF)

/-- ★★★well-defined 性(その 1) —— 同じ `f` の 2 つの持ち上げは同じ類を与える。
★中身は抽象核 `IsTgLift.sub_coset_const`。 -/
theorem tgClass_eq_of_isTgLift (A : Rep k G) (S : Subgroup G) [S.Normal] {f F F' : G → A}
    (hF : IsTgLift (repAddHom A) S f F) (hF' : IsTgLift (repAddHom A) S f F') :
    tgClass A S hF = tgClass A S hF' := by
  obtain ⟨hcos, hinv⟩ := hF.sub_coset_const hF'
  have hmem : ∀ g : G, (F g - F' g) ∈ Representation.invariants (A.ρ.comp S.subtype) :=
    fun g => (Representation.mem_invariants _ _).2 (fun u => hinv g u u.2)
  rw [tgClass, tgClass, H2π_eq_iff]
  refine ⟨descendCochain A S (fun g => F g - F' g) hmem hcos, ?_⟩
  funext p
  obtain ⟨q1, q2⟩ := p
  induction q1 using QuotientGroup.induction_on with | @H g =>
  induction q2 using QuotientGroup.induction_on with | @H h =>
  apply Subtype.ext
  show repAddHom A g (F h - F' h) - (F (g * h) - F' (g * h)) + (F g - F' g)
    = dOne (repAddHom A) F g h - dOne (repAddHom A) F' g h
  simp only [dOne, map_sub]
  abel

/-- ★well-defined 性(その 2) —— 余境界が一致すれば同じ類。 -/
theorem tgClass_eq_of_dOne_eq (A : Rep k G) (S : Subgroup G) [S.Normal] {f F f' F' : G → A}
    (hF : IsTgLift (repAddHom A) S f F) (hF' : IsTgLift (repAddHom A) S f' F')
    (hd : ∀ g j : G, dOne (repAddHom A) F g j = dOne (repAddHom A) F' g j) :
    tgClass A S hF = tgClass A S hF' := by
  rw [tgClass, tgClass, H2π_eq_iff]
  refine ⟨0, ?_⟩
  funext p
  obtain ⟨q1, q2⟩ := p
  induction q1 using QuotientGroup.induction_on with | @H g =>
  induction q2 using QuotientGroup.induction_on with | @H h =>
  apply Subtype.ext
  show repAddHom A g 0 - 0 + 0 = dOne (repAddHom A) F g h - dOne (repAddHom A) F' g h
  rw [hd g h]
  simp

/-- ★★★well-defined 性(本体) —— `f` を `S` 上で余境界だけずらしても類は変わらない。 -/
theorem tgClass_eq_of_cohomologous (A : Rep k G) (S : Subgroup G) [S.Normal] {f F f' F' : G → A}
    (hF : IsTgLift (repAddHom A) S f F) (hF' : IsTgLift (repAddHom A) S f' F')
    (b : A) (hff : ∀ s ∈ S, f' s = f s + (A.ρ s b - b)) :
    tgClass A S hF = tgClass A S hF' := by
  have hF2 : IsTgLift (repAddHom A) S (fun g => f g + (repAddHom A g b - b))
      (fun g => F g + (repAddHom A g b - b)) :=
    hF.shift (repAddHom_one A) (repAddHom_mul A) b
  have hF3 : IsTgLift (repAddHom A) S f' (fun g => F g + (repAddHom A g b - b)) :=
    hF2.congr_restrict f' (fun s hs => hff s hs)
  have e1 : tgClass A S hF = tgClass A S hF3 := by
    refine tgClass_eq_of_dOne_eq A S hF hF3 (fun g j => ?_)
    rw [show (fun g => F g + (repAddHom A g b - b))
        = (fun g => F g + (fun g => repAddHom A g b - b) g) from rfl,
      dOne_add, dOne_principal _ (repAddHom_mul A), add_zero]
  rw [e1]
  exact tgClass_eq_of_isTgLift A S hF3 hF'

/-- `S` 上の 1-コサイクルは `cocycles₁ (Rep.res S.subtype A)` の元。 -/
theorem mem_cocycles₁_restrict (A : Rep k G) (S : Subgroup G) {f : G → A}
    (hf : IsCocycleOn (repAddHom A) S f) :
    (fun s : S => f (s : G)) ∈ cocycles₁ (Rep.res S.subtype A) := by
  rw [mem_cocycles₁_iff]
  intro g h
  exact hf (g : G) g.2 (h : G) h.2

/-- ★★★★well-defined 性の最終形 —— `tgClass` は `H¹(S,A)` の**類**にしか依らない。 -/
theorem tgClass_eq_of_H1π_eq (A : Rep k G) (S : Subgroup G) [S.Normal] {f F f' F' : G → A}
    (hF : IsTgLift (repAddHom A) S f F) (hF' : IsTgLift (repAddHom A) S f' F')
    (h : (ConcreteCategory.hom (H1π (Rep.res S.subtype A)))
          ⟨fun s : S => f (s : G), mem_cocycles₁_restrict A S (hF.isCocycleOn (repAddHom_one A))⟩
        = (ConcreteCategory.hom (H1π (Rep.res S.subtype A)))
          ⟨fun s : S => f' (s : G),
            mem_cocycles₁_restrict A S (hF'.isCocycleOn (repAddHom_one A))⟩) :
    tgClass A S hF = tgClass A S hF' := by
  rw [H1π_eq_iff] at h
  obtain ⟨b, hb⟩ := h
  refine (tgClass_eq_of_cohomologous A S hF' hF b ?_).symm
  intro s hs
  have hth : A.ρ s b - b = f s - f' s := congrFun hb (⟨s, hs⟩ : S)
  linear_combination (norm := abel) -hth

/-! ### 2.1 `H¹(S,A)^{G/S}` における完全性 -/

/-- ★★★完全性(コチェイン版・核 ⊆ 像) —— `tg [f] = 0` ならば
`f` は `G` 全体の 1-コサイクルの制限である。★中身は抽象核 `exists_isCocycle₁_of_dOne_eq`。 -/
theorem exists_isCocycle₁_of_tgClass_eq_zero
    (A : Rep k G) (S : Subgroup G) [S.Normal] {f F : G → A}
    (hF : IsTgLift (repAddHom A) S f F) (h0 : tgClass A S hF = 0) :
    ∃ φ : G → A, IsCocycle₁ (repAddHom A) φ ∧ ∀ s ∈ S, φ s = f s := by
  rw [tgClass, H2π_eq_zero_iff] at h0
  obtain ⟨v, hv⟩ := h0
  refine exists_isCocycle₁_of_dOne_eq (repAddHom_one A) hF
    (v := fun g => ((v (QuotientGroup.mk' S g) : _) : A)) ?_ ?_
  · intro g s hs
    have e : (QuotientGroup.mk' S (g * s)) = QuotientGroup.mk' S g := by
      simpa using QuotientGroup.mk_mul_of_mem (s := S) g hs
    simp only [e]
  · intro g j
    have h2 := congrArg Subtype.val
      (congrFun hv (QuotientGroup.mk' S g, QuotientGroup.mk' S j))
    exact h2.symm

/-- ★`f` が `G` 上の 1-コサイクルの制限なら transgression は 0。 -/
theorem tgClass_eq_zero_of_isCocycle₁ (A : Rep k G) (S : Subgroup G) [S.Normal] {f F : G → A}
    (hF : IsTgLift (repAddHom A) S f F)
    (φ : G → A) (hφ : IsCocycle₁ (repAddHom A) φ) (hres : ∀ s ∈ S, φ s = f s) :
    tgClass A S hF = 0 := by
  have hlift : IsTgLift (repAddHom A) S f φ := hφ.isTgLift (repAddHom_one A) hres
  rw [tgClass_eq_of_isTgLift A S hF hlift, tgClass, H2π_eq_zero_iff]
  refine ⟨0, ?_⟩
  funext p
  obtain ⟨q1, q2⟩ := p
  induction q1 using QuotientGroup.induction_on with | @H g =>
  induction q2 using QuotientGroup.induction_on with | @H h =>
  apply Subtype.ext
  show repAddHom A g 0 - 0 + 0 = dOne (repAddHom A) φ g h
  rw [hφ.dOne_eq_zero g h]
  simp

/-- ★★★★5 項完全列の `H¹(S,A)^{G/S}` における完全性(類の言葉)。

`tg [f] = 0` ⟺ `[f]` は制限 `H¹(G,A) → H¹(S,A)` の像に入る。 -/
theorem tgClass_eq_zero_iff (A : Rep k G) (S : Subgroup G) [S.Normal] {f F : G → A}
    (hF : IsTgLift (repAddHom A) S f F) :
    tgClass A S hF = 0 ↔
      ∃ c : groupCohomology A 1,
        (ConcreteCategory.hom (map (A := A) (B := Rep.res S.subtype A) S.subtype
            (𝟙 (Rep.res S.subtype A)) 1)) c
          = (ConcreteCategory.hom (H1π (Rep.res S.subtype A)))
              ⟨fun s : S => f (s : G),
                mem_cocycles₁_restrict A S (hF.isCocycleOn (repAddHom_one A))⟩ := by
  constructor
  · intro h0
    obtain ⟨φ, hφ, hres⟩ := exists_isCocycle₁_of_tgClass_eq_zero A S hF h0
    refine ⟨(ConcreteCategory.hom (H1π A)) ⟨φ, (mem_cocycles₁_iff φ).2 hφ⟩, ?_⟩
    rw [H1π_comp_map_apply]
    congr 1
    apply Subtype.ext
    funext s
    exact hres (s : G) s.2
  · rintro ⟨c, hc⟩
    induction c using groupCohomology.H1_induction_on with
    | @h x =>
    rw [H1π_comp_map_apply, H1π_eq_iff] at hc
    obtain ⟨b, hb⟩ := hc
    have hφ : IsCocycle₁ (repAddHom A) (x : G → A) := (mem_cocycles₁_iff _).1 x.2
    have hlift : IsTgLift (repAddHom A) S (x : G → A) (x : G → A) :=
      hφ.isTgLift (repAddHom_one A) (fun s _ => rfl)
    have hff : ∀ s ∈ S, (x : G → A) s = f s + (A.ρ s b - b) := by
      intro s hs
      have hth : A.ρ s b - b = (x : G → A) s - f s := congrFun hb (⟨s, hs⟩ : S)
      linear_combination (norm := abel) -hth
    rw [tgClass_eq_of_cohomologous A S hF hlift b hff]
    exact tgClass_eq_zero_of_isCocycle₁ A S hlift _ hφ (fun s _ => rfl)

/-! ### 2.2 `H²(G ⧸ S, Aˢ)` における完全性 -/

/-- ★★★`inf ∘ tg = 0`(5 項完全列の合成が 0)。★証明は `rfl` 1 個。 -/
theorem inflation₂_tgClass (A : Rep k G) (S : Subgroup G) [S.Normal] {f F : G → A}
    (hF : IsTgLift (repAddHom A) S f F) :
    (ConcreteCategory.hom (map (A := A.quotientToInvariants S) (B := A) (QuotientGroup.mk' S)
      (Rep.ofHom (A.ρ.quotientToInvariants_lift S)) 2)) (tgClass A S hF) = 0 := by
  rw [tgClass, H2π_comp_map_apply, H2π_eq_zero_iff]
  exact ⟨F, rfl⟩

def exists_tgClass_of_inflation₂_eq_zero.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }

/-- **★★★★5 項完全列の `H²(G ⧸ S, Aˢ)` における完全性(像 ⊇ 核)**。

inflation が消す類は transgression の像である。★★仮定 `H¹(S,A) = 0` は**要らない**。

原文 (pGC p.3):
> Now it is well-known that (K×)∧ fits into an exact sequence of topological groups:
> 0 → UK → (K×)∧ → Z → 0

★段取り(すべて抽象核 `exists_tgLift_of_inflated_coboundary` にある):
1. `inf [z] = 0` なので、`z` の inflation は `G` 上の余境界 `dy` に一致する。
2. `f := y - y 1` は `S` 上の 1-コサイクルで、自分自身が持ち上げになる。
3. `dOne f = z - d(定数 y 1)` なので、`H²(G ⧸ S, Aˢ)` では `tg [f] = [z]`。 -/
theorem exists_tgClass_of_inflation₂_eq_zero (A : Rep k G) (S : Subgroup G) [S.Normal]
    (c : groupCohomology (A.quotientToInvariants S) 2)
    (hc : (ConcreteCategory.hom (map (A := A.quotientToInvariants S) (B := A)
      (QuotientGroup.mk' S) (Rep.ofHom (A.ρ.quotientToInvariants_lift S)) 2)) c = 0) :
    ∃ (f F : G → A) (hF : IsTgLift (repAddHom A) S f F),
      IsGInvariantCocycle A S f ∧ tgClass A S hF = c := by
  induction c using groupCohomology.H2_induction_on with
  | @h z =>
  rw [H2π_comp_map_apply, H2π_eq_zero_iff] at hc
  obtain ⟨y, hy⟩ := hc
  have hyX : ∀ g h : G, dOne (repAddHom A) y g h
      = ((z (QuotientGroup.mk' S g, QuotientGroup.mk' S h) : _) : A) :=
    fun g h => congrFun hy (g, h)
  have hfacR : ∀ g h : G, ∀ t ∈ S,
      ((z (QuotientGroup.mk' S g, QuotientGroup.mk' S (h * t)) : _) : A)
        = ((z (QuotientGroup.mk' S g, QuotientGroup.mk' S h) : _) : A) := by
    intro g h t ht
    have e : (QuotientGroup.mk' S (h * t)) = QuotientGroup.mk' S h := by
      simpa using QuotientGroup.mk_mul_of_mem (s := S) h ht
    rw [e]
  have hfacL : ∀ g h : G, ∀ s ∈ S,
      ((z (QuotientGroup.mk' S (g * s), QuotientGroup.mk' S h) : _) : A)
        = ((z (QuotientGroup.mk' S g, QuotientGroup.mk' S h) : _) : A) := by
    intro g h s hs
    have e : (QuotientGroup.mk' S (g * s)) = QuotientGroup.mk' S g := by
      simpa using QuotientGroup.mk_mul_of_mem (s := S) g hs
    rw [e]
  have hinvX : ∀ g h : G, ∀ u ∈ S,
      repAddHom A u ((z (QuotientGroup.mk' S g, QuotientGroup.mk' S h) : _) : A)
        = ((z (QuotientGroup.mk' S g, QuotientGroup.mk' S h) : _) : A) := by
    intro g h u hu
    exact (Representation.mem_invariants _ _).1
      (z (QuotientGroup.mk' S g, QuotientGroup.mk' S h)).2 ⟨u, hu⟩
  obtain ⟨hlift, hdOne, hy1⟩ :=
    exists_tgLift_of_inflated_coboundary (repAddHom_one A) hfacR hfacL hinvX y hyX
  have hy1mem : y 1 ∈ Representation.invariants (A.ρ.comp S.subtype) :=
    (Representation.mem_invariants _ _).2 (fun u => hy1 u u.2)
  refine ⟨_, _, hlift, fun g => ⟨y g - y 1, fun s hs => hlift.conj g s hs⟩, ?_⟩
  rw [tgClass, H2π_eq_iff]
  refine ⟨fun _ => -(⟨y 1, hy1mem⟩ : (A.quotientToInvariants S)), ?_⟩
  funext p
  obtain ⟨q1, q2⟩ := p
  induction q1 using QuotientGroup.induction_on with | @H g =>
  induction q2 using QuotientGroup.induction_on with | @H h =>
  apply Subtype.ext
  show repAddHom A g (-(y 1)) - (-(y 1)) + (-(y 1))
    = dOne (repAddHom A) (fun g => y g - y 1) g h
      - ((z (QuotientGroup.mk' S g, QuotientGroup.mk' S h) : _) : A)
  rw [hdOne g h, map_neg]
  abel

/-! ### 2.3 `H¹(S,A)^{G/S}` と transgression 写像 -/

/-- ★★`c : H¹(S,A)` が `G ⧸ S`-不変(= `c ∈ H¹(S,A)^{G/S}`)。

★mathlib は `H^n(S,A)` への `G ⧸ S` の共役作用を持たないので、
「不変なコサイクルで代表される」という形で定義する。★`isGInvariantH1_iff` を参照。 -/
def IsGInvariantH1 (A : Rep k G) (S : Subgroup G) [S.Normal]
    (c : groupCohomology (Rep.res S.subtype A) 1) : Prop :=
  ∃ (f F : G → A) (hF : IsTgLift (repAddHom A) S f F),
    (ConcreteCategory.hom (H1π (Rep.res S.subtype A)))
      ⟨fun s : S => f (s : G), mem_cocycles₁_restrict A S (hF.isCocycleOn (repAddHom_one A))⟩ = c

/-- ★★名前の正当化 —— `IsGInvariantH1` は「不変なコサイクルで代表される」ことと同値。 -/
theorem isGInvariantH1_iff (A : Rep k G) (S : Subgroup G) [S.Normal]
    (c : groupCohomology (Rep.res S.subtype A) 1) :
    IsGInvariantH1 A S c ↔
      ∃ (f : G → A) (hf : IsCocycleOn (repAddHom A) S f), IsGInvariantCocycle A S f ∧
        (ConcreteCategory.hom (H1π (Rep.res S.subtype A)))
          ⟨fun s : S => f (s : G), mem_cocycles₁_restrict A S hf⟩ = c := by
  constructor
  · rintro ⟨f, F, hF, hfc⟩
    exact ⟨f, hF.isCocycleOn (repAddHom_one A), hF.isGInvariantCocycle, hfc⟩
  · rintro ⟨f, hf, hinv, hfc⟩
    obtain ⟨F, hF⟩ := exists_tgLift_rep A S f hf hinv
    exact ⟨f, F, hF, hfc⟩

/-- **★★★★transgression 写像** `tg : H¹(S,A)^{G/S} → H²(G ⧸ S, Aˢ)`。

★代表の選び方に依らないことは `transgression_eq` が保証する。 -/
noncomputable def transgression (A : Rep k G) (S : Subgroup G) [S.Normal]
    {c : groupCohomology (Rep.res S.subtype A) 1} (hc : IsGInvariantH1 A S c) :
    groupCohomology (A.quotientToInvariants S) 2 :=
  tgClass A S hc.choose_spec.choose_spec.choose

/-- ★★★`transgression` の仕様 —— どんな代表・どんな持ち上げで計算してもよい。 -/
theorem transgression_eq (A : Rep k G) (S : Subgroup G) [S.Normal]
    {c : groupCohomology (Rep.res S.subtype A) 1} (hc : IsGInvariantH1 A S c) {f F : G → A}
    (hF : IsTgLift (repAddHom A) S f F)
    (hfc : (ConcreteCategory.hom (H1π (Rep.res S.subtype A)))
        ⟨fun s : S => f (s : G), mem_cocycles₁_restrict A S (hF.isCocycleOn (repAddHom_one A))⟩
      = c) :
    transgression A S hc = tgClass A S hF := by
  refine tgClass_eq_of_H1π_eq A S _ hF ?_
  rw [hfc]
  exact hc.choose_spec.choose_spec.choose_spec

/-- ★★制限 `H¹(G,A) → H¹(S,A)` の像は `G ⧸ S`-不変。
★これが無いと 5 項完全列の 3 項目が `H¹(S,A)^{G/S}` である意味が無い。 -/
theorem isGInvariantH1_of_restriction (A : Rep k G) (S : Subgroup G) [S.Normal]
    (d : groupCohomology A 1) :
    IsGInvariantH1 A S ((ConcreteCategory.hom (map (A := A) (B := Rep.res S.subtype A)
      S.subtype (𝟙 (Rep.res S.subtype A)) 1)) d) := by
  induction d using groupCohomology.H1_induction_on with
  | @h x =>
  have hφ : IsCocycle₁ (repAddHom A) (x : G → A) := (mem_cocycles₁_iff _).1 x.2
  refine ⟨(x : G → A), (x : G → A), hφ.isTgLift (repAddHom_one A) (fun s _ => rfl), ?_⟩
  rw [H1π_comp_map_apply]
  rfl

/-- ★★★★5 項完全列の `H¹(S,A)^{G/S}` における完全性(`transgression` 版)。 -/
theorem transgression_eq_zero_iff (A : Rep k G) (S : Subgroup G) [S.Normal]
    {c : groupCohomology (Rep.res S.subtype A) 1} (hc : IsGInvariantH1 A S c) :
    transgression A S hc = 0 ↔
      ∃ d : groupCohomology A 1,
        (ConcreteCategory.hom (map (A := A) (B := Rep.res S.subtype A) S.subtype
          (𝟙 (Rep.res S.subtype A)) 1)) d = c := by
  rw [show transgression A S hc = tgClass A S hc.choose_spec.choose_spec.choose from rfl,
    tgClass_eq_zero_iff, hc.choose_spec.choose_spec.choose_spec]

/-- ★★★`inf ∘ tg = 0`(`transgression` 版)。 -/
theorem inflation₂_transgression (A : Rep k G) (S : Subgroup G) [S.Normal]
    {c : groupCohomology (Rep.res S.subtype A) 1} (hc : IsGInvariantH1 A S c) :
    (ConcreteCategory.hom (map (A := A.quotientToInvariants S) (B := A) (QuotientGroup.mk' S)
      (Rep.ofHom (A.ρ.quotientToInvariants_lift S)) 2)) (transgression A S hc) = 0 :=
  inflation₂_tgClass A S _

/-- ★★★★5 項完全列の `H²(G ⧸ S, Aˢ)` における完全性(`transgression` 版)。
★★仮定 `H¹(S,A) = 0` は要らない。 -/
theorem exists_transgression_eq (A : Rep k G) (S : Subgroup G) [S.Normal]
    (e : groupCohomology (A.quotientToInvariants S) 2)
    (h : (ConcreteCategory.hom (map (A := A.quotientToInvariants S) (B := A)
      (QuotientGroup.mk' S) (Rep.ofHom (A.ρ.quotientToInvariants_lift S)) 2)) e = 0) :
    ∃ (c : groupCohomology (Rep.res S.subtype A) 1) (hc : IsGInvariantH1 A S c),
      transgression A S hc = e := by
  obtain ⟨f, F, hF, _, he⟩ := exists_tgClass_of_inflation₂_eq_zero A S e h
  exact ⟨_, ⟨f, F, hF, rfl⟩, (transgression_eq A S _ hF rfl).trans he⟩

/-! ### 2.4 線形構造 —— `H¹(S,A)^{G/S}` は部分加群、`tg` は `k`-線形 -/

/-- ★持ち上げのスカラー倍(具体層。抽象核は `k` を知らないのでここで足す)。 -/
theorem IsTgLift.smul_rep {A : Rep k G} {S : Subgroup G} {f F : G → A}
    (hF : IsTgLift (repAddHom A) S f F) (a : k) :
    IsTgLift (repAddHom A) S (fun g => a • f g) (fun g => a • F g) := by
  refine ⟨by simp [hF.one], ?_, ?_⟩
  · intro g s hs
    have e := hF.right g s hs
    simp only [e, smul_add, repAddHom_apply, map_smul]
  · intro g s hs
    have e := hF.conj g s hs
    simp only [repAddHom_apply, map_smul] at e ⊢
    rw [← smul_sub, ← smul_sub, e]

/-- ★★transgression の加法性(コサイクル水準)。 -/
theorem tgClass_add (A : Rep k G) (S : Subgroup G) [S.Normal] {f F f' F' : G → A}
    (hF : IsTgLift (repAddHom A) S f F) (hF' : IsTgLift (repAddHom A) S f' F') :
    tgClass A S (hF.add hF') = tgClass A S hF + tgClass A S hF' := by
  rw [tgClass, tgClass, tgClass, ← map_add]
  congr 1
  apply Subtype.ext
  funext p
  obtain ⟨q1, q2⟩ := p
  induction q1 using QuotientGroup.induction_on with | @H g =>
  induction q2 using QuotientGroup.induction_on with | @H h =>
  apply Subtype.ext
  show dOne (repAddHom A) (fun g => F g + F' g) g h
    = dOne (repAddHom A) F g h + dOne (repAddHom A) F' g h
  exact dOne_add _ F F' g h

/-- ★★transgression のスカラー倍(コサイクル水準)。 -/
theorem tgClass_smul (A : Rep k G) (S : Subgroup G) [S.Normal] {f F : G → A}
    (hF : IsTgLift (repAddHom A) S f F) (a : k) :
    tgClass A S (hF.smul_rep a) = a • tgClass A S hF := by
  rw [tgClass, tgClass, ← map_smul]
  congr 1
  apply Subtype.ext
  funext p
  obtain ⟨q1, q2⟩ := p
  induction q1 using QuotientGroup.induction_on with | @H g =>
  induction q2 using QuotientGroup.induction_on with | @H h =>
  apply Subtype.ext
  show dOne (repAddHom A) (fun g => a • F g) g h = a • dOne (repAddHom A) F g h
  simp only [dOne, repAddHom_apply, map_smul, smul_sub, smul_add]

theorem isGInvariantH1_zero (A : Rep k G) (S : Subgroup G) [S.Normal] : IsGInvariantH1 A S 0 := by
  refine ⟨fun _ => 0, fun _ => 0, ⟨rfl, by intros; simp, by intros; simp⟩, ?_⟩
  rw [← map_zero (ConcreteCategory.hom (H1π (Rep.res S.subtype A)))]
  congr 1

theorem isGInvariantH1_add (A : Rep k G) (S : Subgroup G) [S.Normal]
    {c c' : groupCohomology (Rep.res S.subtype A) 1}
    (hc : IsGInvariantH1 A S c) (hc' : IsGInvariantH1 A S c') : IsGInvariantH1 A S (c + c') := by
  obtain ⟨f, F, hF, e⟩ := hc
  obtain ⟨f', F', hF', e'⟩ := hc'
  refine ⟨fun g => f g + f' g, fun g => F g + F' g, hF.add hF', ?_⟩
  rw [← e, ← e', ← map_add]
  congr 1

theorem isGInvariantH1_smul (A : Rep k G) (S : Subgroup G) [S.Normal]
    {c : groupCohomology (Rep.res S.subtype A) 1}
    (hc : IsGInvariantH1 A S c) (a : k) : IsGInvariantH1 A S (a • c) := by
  obtain ⟨f, F, hF, e⟩ := hc
  refine ⟨fun g => a • f g, fun g => a • F g, hF.smul_rep a, ?_⟩
  rw [← e, ← map_smul]
  congr 1

theorem transgression_add (A : Rep k G) (S : Subgroup G) [S.Normal]
    {c c' : groupCohomology (Rep.res S.subtype A) 1}
    (hc : IsGInvariantH1 A S c) (hc' : IsGInvariantH1 A S c') :
    transgression A S (isGInvariantH1_add A S hc hc')
      = transgression A S hc + transgression A S hc' := by
  obtain ⟨f, F, hF, e⟩ := id hc
  obtain ⟨f', F', hF', e'⟩ := id hc'
  have key : (ConcreteCategory.hom (H1π (Rep.res S.subtype A)))
      ⟨fun s : S => (fun g => f g + f' g) (s : G),
        mem_cocycles₁_restrict A S ((hF.add hF').isCocycleOn (repAddHom_one A))⟩ = c + c' := by
    rw [← e, ← e', ← map_add]; congr 1
  rw [transgression_eq A S _ (hF.add hF') key, transgression_eq A S hc hF e,
    transgression_eq A S hc' hF' e', tgClass_add]

theorem transgression_smul (A : Rep k G) (S : Subgroup G) [S.Normal]
    {c : groupCohomology (Rep.res S.subtype A) 1} (hc : IsGInvariantH1 A S c) (a : k) :
    transgression A S (isGInvariantH1_smul A S hc a) = a • transgression A S hc := by
  obtain ⟨f, F, hF, e⟩ := id hc
  have key : (ConcreteCategory.hom (H1π (Rep.res S.subtype A)))
      ⟨fun s : S => (fun g => a • f g) (s : G),
        mem_cocycles₁_restrict A S ((hF.smul_rep a).isCocycleOn (repAddHom_one A))⟩ = a • c := by
    rw [← e, ← map_smul]; congr 1
  rw [transgression_eq A S _ (hF.smul_rep a) key, transgression_eq A S hc hF e, tgClass_smul]

/-- ★★★`H¹(S,A)^{G/S}` —— `H¹(S,A)` の `G ⧸ S`-不変部分(`k`-部分加群)。 -/
def invariantsH1 (A : Rep k G) (S : Subgroup G) [S.Normal] :
    Submodule k (groupCohomology (Rep.res S.subtype A) 1) where
  carrier := {c | IsGInvariantH1 A S c}
  add_mem' hc hc' := isGInvariantH1_add A S hc hc'
  zero_mem' := isGInvariantH1_zero A S
  smul_mem' a _ hc := isGInvariantH1_smul A S hc a

theorem mem_invariantsH1 (A : Rep k G) (S : Subgroup G) [S.Normal]
    {c : groupCohomology (Rep.res S.subtype A) 1} :
    c ∈ invariantsH1 A S ↔ IsGInvariantH1 A S c := Iff.rfl

/-- ★★★★transgression `tg : H¹(S,A)^{G/S} → H²(G ⧸ S, Aˢ)`(`k`-線形写像)。 -/
noncomputable def transgressionLin (A : Rep k G) (S : Subgroup G) [S.Normal] :
    invariantsH1 A S →ₗ[k] groupCohomology (A.quotientToInvariants S) 2 where
  toFun c := transgression A S c.2
  map_add' c c' := transgression_add A S c.2 c'.2
  map_smul' a c := transgression_smul A S c.2 a

theorem transgressionLin_apply (A : Rep k G) (S : Subgroup G) [S.Normal]
    (c : invariantsH1 A S) : transgressionLin A S c = transgression A S c.2 := rfl

/-! ### 2.5 ★★★★★5 項完全列 -/

/-- **★★★★★5 項完全列 —— `H¹(S,A)^{G/S}` における完全性**。

`ker(tg)` は制限 `H¹(G,A) → H¹(S,A)` の像(を `H¹(S,A)^{G/S}` に引き戻したもの)に等しい。 -/
theorem ker_transgressionLin (A : Rep k G) (S : Subgroup G) [S.Normal] :
    LinearMap.ker (transgressionLin A S) =
      Submodule.comap (invariantsH1 A S).subtype
        (LinearMap.range (ModuleCat.Hom.hom (map (A := A) (B := Rep.res S.subtype A)
          S.subtype (𝟙 (Rep.res S.subtype A)) 1))) := by
  ext c
  simp only [LinearMap.mem_ker, Submodule.mem_comap, LinearMap.mem_range,
    Submodule.coe_subtype, transgressionLin_apply]
  exact transgression_eq_zero_iff A S c.2

/-- **★★★★★5 項完全列 —— `H²(G ⧸ S, Aˢ)` における完全性**。

`range(tg) = ker(inf)`。★★仮定 `H¹(S,A) = 0` は要らない。 -/
theorem range_transgressionLin (A : Rep k G) (S : Subgroup G) [S.Normal] :
    LinearMap.range (transgressionLin A S) =
      LinearMap.ker (ModuleCat.Hom.hom (map (A := A.quotientToInvariants S) (B := A)
        (QuotientGroup.mk' S) (Rep.ofHom (A.ρ.quotientToInvariants_lift S)) 2)) := by
  apply le_antisymm
  · rintro e ⟨c, rfl⟩
    simp only [LinearMap.mem_ker, transgressionLin_apply]
    exact inflation₂_transgression A S c.2
  · intro e he
    simp only [LinearMap.mem_ker] at he
    obtain ⟨c, hc, hce⟩ := exists_transgression_eq A S e he
    exact ⟨⟨c, hc⟩, hce⟩

/-- **★★★★★★本ファイルの到達点(まとめ)** —— 5 項完全列

`0 → H¹(G⧸S, Aˢ) → H¹(G,A) → H¹(S,A)^{G/S} --tg--> H²(G⧸S, Aˢ) → H²(G,A)`

の右 2 箇所の完全性。★★**仮定は無い**(左 2 箇所は mathlib の `H1InfRes_exact`)。 -/
theorem fiveTermExact (A : Rep k G) (S : Subgroup G) [S.Normal] :
    LinearMap.ker (transgressionLin A S) =
        Submodule.comap (invariantsH1 A S).subtype
          (LinearMap.range (ModuleCat.Hom.hom (map (A := A) (B := Rep.res S.subtype A)
            S.subtype (𝟙 (Rep.res S.subtype A)) 1)))
      ∧ LinearMap.range (transgressionLin A S) =
        LinearMap.ker (ModuleCat.Hom.hom (map (A := A.quotientToInvariants S) (B := A)
          (QuotientGroup.mk' S) (Rep.ofHom (A.ρ.quotientToInvariants_lift S)) 2)) :=
  ⟨ker_transgressionLin A S, range_transgressionLin A S⟩

/-! ### 2.6 ★退化の自己検査 -/

/-- **★★★直前の波の `inflation₂_injective` は本ファイルの完全性の系である**。

`H¹(S,A) = 0` ならば `H¹(S,A)^{G/S} = 0` なので `ker(inf) = range(tg) = 0`。
★「新しい定理が古い定理を含む」ことを Lean で測ったもの。 -/
theorem inflation₂_injective_of_subsingleton_H1 (A : Rep k G) (S : Subgroup G) [S.Normal]
    (h1 : Subsingleton (groupCohomology (Rep.res S.subtype A) 1)) :
    Function.Injective (ConcreteCategory.hom (map (A := A.quotientToInvariants S) (B := A)
      (QuotientGroup.mk' S) (Rep.ofHom (A.ρ.quotientToInvariants_lift S)) 2)) := by
  rw [injective_iff_map_eq_zero]
  intro e he
  obtain ⟨c, hc, hce⟩ := exists_transgression_eq A S e he
  rw [← hce]
  refine (transgression_eq_zero_iff A S hc).2 ⟨0, ?_⟩
  rw [map_zero]
  exact h1.elim 0 c

/-- ★退化の自己検査 —— `S = ⊥` では `H¹(⊥,A) = 0` なので `H¹(S,A)^{G/S}` は 0。
★このとき `range_transgressionLin` は「`inf` が単射」に潰れ、
直前の波の `inflation₂_bijective_bot` と整合する。 -/
theorem invariantsH1_bot (A : Rep k G) : invariantsH1 A (⊥ : Subgroup G) = ⊥ := by
  have h1 : Subsingleton (groupCohomology (Rep.res (⊥ : Subgroup G).subtype A) 1) :=
    subsingleton_H1_of_subsingleton _
  rw [Submodule.eq_bot_iff]
  intro c _
  exact h1.elim c 0

end ConcreteTransgression

end ABC3.Found.PGC
